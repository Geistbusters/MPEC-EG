c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdsuc1(prod,ncsfs,ncsfd,ncsf2,xint,vs,gs,vd,gd,lsa,   6d10s21
     $     jsbv,isbv12,nvirt,multh,nrootu,nsymb,sr2,mrowx,ldebug,loop,
     $     loopx,bc,ibc)                                                11d14s22
      implicit real*8 (a-h,o-z)
      logical ldebug
      dimension prod(ncsfs,*),ncsf2(*),xint(*),vs(*),gs(*),vd(*),gd(*), 6d10s21
     $     nvirt(*),multh(8,8),nfdat(2),intden(2),itrans3(2),itrans(2)  6d11s21
      include "common.store"                                            6d10s21
      nfdat(1)=ncsf2(1)                                                 6d10s21
      nfdat(2)=ncsfd-ncsf2(1)                                           6d10s21
      iaddd=1                                                           6d11s21
      do l=1,2                                                          6d10s21
       nrow=ncsfs*nfdat(l)                                              6d10s21
       itrans(l)=ibcoff                                                 6d10s21
       ibcoff=itrans(l)+nvirt(lsa)*nrow                                     6d10s21
       call enough(' hcdsuc1.  1',bc,ibc)
       do id=0,nfdat(l)-1                                               6d11s21
        idp=id+iaddd                                                    6d11s21
        do is=1,ncsfs                                                   6d11s21
         jtrans=itrans(l)-1+nvirt(lsa)*(is-1+ncsfs*id)                  6d11s21
         do jv=1,nvirt(lsa)                                               6d10s21
          bc(jtrans+jv)=xint(jv)*prod(is,idp)                           6d11s21
         end do                                                         6d11s21
        end do                                                           6d10s21
       end do                                                            6d10s21
       iaddd=iaddd+nfdat(l)                                             6d11s21
      end do                                                            6d11s21
      if(isbv12.eq.1)then                                               6d10s21
       nrow=ncsfs*nfdat(1)                                              6d11s21
       itrans2=ibcoff                                                   6d10s21
       ibcoff=itrans2+nrow*nvirt(lsa)                                   6d10s21
       call enough(' hcdsuc1.  2',bc,ibc)
       do i=0,nrow*nvirt(lsa)-1                                         6d10s21
        bc(itrans2+i)=bc(itrans(1)+i)*sr2                               6d11s21
       end do                                                           6d10s21
      else                                                              6d10s21
       ivst=ibcoff                                                      6d10s21
       ibcoff=ivst+nvirt(jsbv)*ncsfs*nrootu                             6d10s21
       call enough(' hcdsuc1.  3',bc,ibc)
       do is=0,ncsfs-1                                                  6d10s21
        do ir=0,nrootu-1                                                6d10s21
         iadf=1+nvirt(jsbv)*(ir+nrootu*is)                              6d10s21
         do iv=0,nvirt(jsbv)-1                                           6d10s21
          iadt=ivst+is+ncsfs*(iv+nvirt(jsbv)*ir)                        6d10s21
          bc(iadt)=vs(iadf+iv)                                           6d10s21
         end do                                                         6d11s21
        end do                                                          6d10s21
       end do                                                           6d10s21
       iaddd=0                                                          6d10s21
       do l=1,2                                                         6d10s21
        intden(l)=ibcoff                                                6d10s21
        itrans3(l)=intden(l)+nvirt(lsa)*ncsfs*nfdat(l)                  6d10s21
        ibcoff=itrans3(l)+nfdat(l)*nvirt(lsa)*ncsfs                     6d10s21
        call enough(' hcdsuc1.  4',bc,ibc)
        do id=0,nfdat(l)-1                                              6d10s21
         idp=id+iaddd                                                   6d10s21
         do is=0,ncsfs-1                                                6d10s21
          do iv=0,nvirt(lsa)-1                                          6d10s21
           iad1=itrans(l)+iv+nvirt(lsa)*(is+ncsfs*id)                   6d11s21
           iad2=intden(l)+is+ncsfs*(id+nfdat(l)*iv)                     6d10s21
           bc(iad2)=bc(iad1)                                              6d10s21
           iad3=itrans3(l)+id+nfdat(l)*(iv+nvirt(lsa)*is)               6d10s21
           bc(iad3)=bc(iad1)                                            6d10s21
          end do                                                        6d10s21
         end do                                                         6d10s21
        end do                                                          6d10s21
        iaddd=iaddd+nfdat(l)                                            6d10s21
       end do                                                           6d10s21
      end if                                                            6d10s21
      ioffvd=1                                                          6d10s21
      do isbv1=1,nsymb                                                  6d10s21
       isbv2=multh(isbv1,isbv12)                                        6d10s21
       if(isbv2.ge.isbv1)then                                           6d10s21
        if(isbv1.eq.isbv2)then                                          6d10s21
         if(jsbv.eq.isbv1.and.jsbv.eq.lsa)then                          6d10s21
          do id=0,ncsf2(1)-1                                            6d10s21
           do is=0,ncsfs-1                                              6d14s21
            jtrans=itrans2+nvirt(jsbv)*(is+ncsfs*id)                     6d10s21
            jdd=ioffvd+nvirt(jsbv)*nrootu*id                            6d10s21
            jss=1+nvirt(jsbv)*nrootu*is                                 6d10s21
            do ir=0,nrootu-1                                            6d10s21
             do jv=0,nvirt(jsbv)-1                                      6d10s21
              gs(jss+jv)=gs(jss+jv)+bc(jtrans+jv)*vd(jdd+jv)            6d15s21
              gd(jdd+jv)=gd(jdd+jv)+bc(jtrans+jv)*vs(jss+jv)            6d15s21
             end do                                                     6d10s21
             jss=jss+nvirt(jsbv)                                        6d10s21
             jdd=jdd+nvirt(jsbv)                                        6d10s21
            end do                                                      6d10s21
           end do                                                       6d10s21
          end do                                                        6d10s21
         end if                                                         6d10s21
         nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                          6d10s21
         ioffvd=ioffvd+nrootu*nvirt(isbv1)*ncsf2(1)                     6d10s21
         isw=0                                                          6d10s21
        else                                                            6d10s21
         nvv=nvirt(isbv1)*nvirt(isbv2)                                  6d10s21
         isw=1                                                          6d10s21
        end if                                                          6d10s21
        if((jsbv.eq.isbv1.and.isbv2.eq.lsa).or.                         6d10s21
     $     (jsbv.eq.isbv2.and.isbv1.eq.lsa))then                        6d10s21
         if(isbv1.eq.jsbv)then                                          6d10s21
          isbvu=isbv2                                                   6d10s21
          tf=1d0                                                        6d10s21
         else                                                           6d10s21
          isbvu=isbv1                                                   6d10s21
          tf=-1d0                                                       6d10s21
         end if                                                         6d10s21
         do l=1,2                                                       6d10s21
          if(l.eq.1)then                                                6d10s21
           factt=1d0                                                    6d10s21
           tf2=1d0                                                      6d10s21
           iaddd=0                                                      6d10s21
          else                                                          6d10s21
           iaddd=nfdat(1)                                               6d10s21
           factt=tf                                                     6d10s21
           tf2=-1d0                                                     6d10s21
          end if                                                        6d10s21
          if(isbv1.eq.isbv2)then
           do id=0,nfdat(l)-1                                           6d10s21
            idp=id+iaddd                                                6d10s21
            fa=factt                                                    6d10s21
            fb=factt*tf2                                                6d10s21
            do is=0,ncsfs-1                                             6d10s21
             do ir=0,nrootu-1                                           6d10s21
              iadv=ioffvd+nvv*(ir+nrootu*idp)                           6d10s21
              iadt=itrans(l)+nvirt(jsbv)*(is+ncsfs*id)                  6d10s21
              iasv=1+nvirt(jsbv)*(ir+nrootu*is)                         6d10s21
              do iv2=0,nvirt(isbv2)-1                                   6d10s21
               do iv1=0,iv2-1                                           6d10s21
                ff=vd(iadv+iv1)*tf2                                     6d10s21
                gs(iasv+iv2)=gs(iasv+iv2)+bc(iadt+iv1)*vd(iadv+iv1)*fb  6d10s21
                gs(iasv+iv1)=gs(iasv+iv1)+bc(iadt+iv2)*vd(iadv+iv1)*fa  6d10s21
                gd(iadv+iv1)=gd(iadv+iv1)+bc(iadt+iv1)*vs(iasv+iv2)*fb  6d10s21
                gd(iadv+iv1)=gd(iadv+iv1)+bc(iadt+iv2)*vs(iasv+iv1)*fa  6d10s21
               end do                                                   6d10s21
               iadv=iadv+iv2                                            6d10s21
              end do                                                    6d10s21
             end do                                                     6d10s21
            end do                                                      6d10s21
           end do                                                       6d10s21
          else if(lsa.eq.isbv1)then                                     6d10s21
c     r v2 j, j f v1, v1v2 r f
c     f v1 r v2, f v1 j, j r v2
           ivdt=ibcoff                                                  6d10s21
           ibcoff=ivdt+nvv*nrootu*nfdat(l)                              6d10s21
           call enough(' hcdsuc1.  5',bc,ibc)
           do id=0,nfdat(l)-1                                           6d10s21
            idp=id+iaddd                                                6d10s21
            do ir=0,nrootu-1                                            6d10s21
             do iv2=0,nvirt(isbv2)-1                                    6d10s21
              do iv1=0,nvirt(isbv1)-1                                   6d10s21
               iad1=ioffvd+iv1+nvirt(isbv1)*(iv2+nvirt(isbv2)*(ir       6d10s21
     $              +nrootu*idp))                                       6d10s21
               iad2=ivdt+id+nfdat(l)*(iv1+nvirt(isbv1)*(iv2             6d10s21
     $              +nvirt(isbv2)*ir))                                  6d10s21
               bc(iad2)=vd(iad1)                                        6d10s21
              end do                                                    6d10s21
             end do
            end do
           end do                                                       6d10s21
           mcol=nrootu*nvirt(jsbv)                                      6d10s21
           mrow=nvirt(isbv1)*nfdat(l)                                   6d10s21
           itmpsg=ibcoff                                                6d10s21
           if(min(mcol,mrow).gt.0)then                                  8d19s21
            ibcoff=itmpsg+mcol*ncsfs                                     6d10s21
            call enough(' hcdsuc1.  6',bc,ibc)
            call dgemm('n','n',ncsfs,mcol,mrow,factt,bc(intden(l)),     8d16s21
     $                   ncsfs,bc(ivdt),mrow,0d0,bc(itmpsg),ncsfs)      8d16s21
            jtmpsg=itmpsg                                                6d10s21
            do ir=0,nrootu-1                                             6d10s21
             do jv=1,nvirt(jsbv)                                         6d10s21
              do j=0,ncsfs-1                                             6d10s21
               iads=jv+nvirt(jsbv)*(ir+nrootu*j)                         6d10s21
               gs(iads)=gs(iads)+bc(jtmpsg+j)                            6d10s21
              end do                                                     6d10s21
              jtmpsg=jtmpsg+ncsfs                                        6d10s21
             end do                                                      6d10s21
            end do                                                       6d10s21
           end if                                                       8d16s21
           ibcoff=ivdt                                                  6d10s21
           if(min(mcol,mrow).gt.0)then                                  8d19s21
            itmpdg=ibcoff                                                6d10s21
            ibcoff=itmpdg+mrow*mcol                                      6d10s21
            call enough(' hcdsuc1.  7',bc,ibc)
            call dgemm('n','n',mrow,mcol,ncsfs,factt,bc(itrans3(l)),    8d16s21
     $                   mrow,bc(ivst),ncsfs,0d0,bc(itmpdg),mrow)       8d16s21
            jtmpdg=itmpdg                                                6d10s21
            do ir=0,nrootu-1                                             6d10s21
             do iv2=0,nvirt(isbv2)-1                                     6d10s21
              do iv1=0,nvirt(isbv1)-1                                    6d10s21
               do i=0,nfdat(l)-1                                         6d10s21
                iad=ioffvd+iv1+nvirt(isbv1)*(iv2+nvirt(isbv2)*(ir        6d10s21
     $              +nrootu*(i+iaddd)))                                 6d10s21
                gd(iad)=gd(iad)+bc(jtmpdg+i)                                     6d10s21
               end do                                                    6d10s21
               jtmpdg=jtmpdg+nfdat(l)                                    6d10s21
              end do                                                     6d10s21
             end do                                                      6d10s21
            end do                                                       6d10s21
            ibcoff=itmpdg                                                6d10s21
           end if                                                       8d16s21
          else                                                          6d10s21
c     r v1 j, j f v2, v1v2 r f
c     f v1 r v2, f v2 j, j r v1
           mrow=nrootu*nvirt(isbv1)                                     6d10s21
           mcol=nfdat(l)*nvirt(isbv2)                                   6d10s21
           itmpdv=ibcoff                                                6d10s21
           itmpgs=itmpdv+mrow*mcol                                      6d10s21
           ibcoff=itmpgs+ncsfs*mrow                                     6d10s21
           call enough(' hcdsuc1.  8',bc,ibc)
           do k=0,nfdat(l)-1                                            6d10s21
            kk=k+iaddd                                                  6d10s21
            do ir=0,nrootu-1                                            6d10s21
             do iv2=0,nvirt(isbv2)-1                                    6d10s21
              do iv1=0,nvirt(isbv1)-1                                   6d10s21
               iad1=ioffvd+iv1+nvirt(isbv1)*(iv2                        6d10s21
     $                        +nvirt(isbv2)*(ir+nrootu*kk))             6d10s21
               iad2=itmpdv+k+nfdat(l)*(iv2+nvirt(isbv2)                 6d10s21
     $                       *(iv1+nvirt(isbv1)*ir))                    6d10s21
               bc(iad2)=vd(iad1)                                        6d10s21
              end do                                                    6d10s21
             end do                                                     6d10s21
            end do                                                      6d10s21
           end do                                                       6d10s21
           if(min(mrow,mcol).gt.0)then                                  8d19s21
            call dgemm('n','n',ncsfs,mrow,mcol,factt,                    6d10s21
     $          bc(intden(l)),ncsfs,bc(itmpdv),mcol,0d0,                6d10s21
     $          bc(itmpgs),ncsfs)                                       6d10s21
c     intden is,id,iv
c     itmpdv id,iv2,iv1,r
            do i=0,mrow-1                                                6d10s21
             do j=0,ncsfs-1                                              6d10s21
              ji=itmpgs+j+ncsfs*i                                        6d10s21
              ij=i+1+mrow*j                                              6d10s21
              gs(ij)=gs(ij)+bc(ji)                                       6d10s21
             end do                                                      6d10s21
            end do                                                       6d10s21
            call dgemm('n','n',mcol,mrow,ncsfs,factt,                    6d10s21
     $          bc(itrans3(l)),mcol,bc(ivst),ncsfs,0d0,                 6d15s21
     $          bc(itmpdv),mcol)                                        6d10s21
c     itrans iv,is,id
c     intden is,id,iv
c     itrans3 id,iv,is
c     ivst is,jv,r
            do ir=0,nrootu-1                                               6d10s21
             do iv1=0,nvirt(isbv1)-1                                     6d10s21
              do iv2=0,nvirt(isbv2)-1                                    6d10s21
               iad1=itmpdv+nfdat(l)*(iv2+nvirt(isbv2)                    6d10s21
     $                       *(iv1+nvirt(isbv1)*ir))                    6d10s21
               do k=0,nfdat(l)-1                                         6d10s21
                kk=k+iaddd                                               6d10s21
                iad2=ioffvd+iv1+nvirt(isbv1)*(iv2                        6d10s21
     $                       +nvirt(isbv2)*(ir+nrootu*kk))              6d10s21
                gd(iad2)=gd(iad2)+bc(iad1+k)                             6d10s21
               end do                                                    6d10s21
              end do                                                     6d10s21
             end do                                                      6d10s21
            end do                                                       6d10s21
           end if                                                       8d16s21
           ibcoff=itmpdv                                                6d10s21
          end if                                                        6d10s21
         end do                                                         6d10s21
        end if                                                          6d10s21
        ioffvd=ioffvd+nrootu*nvv*ncsfd                                  6d10s21
       end if                                                           6d10s21
      end do                                                            6d10s21
      return                                                            6d10s21
      end                                                               6d10s21
