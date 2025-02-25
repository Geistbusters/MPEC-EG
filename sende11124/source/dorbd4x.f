c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dorbd4x(nsymb,idoubo,irefo,noc,nvirt,nsdlk1,isblk1,i3x,11d17s23
     $     itt,nbasdws,isymmrci,multh,vdinout,nfdat,srh,ndoub,nroot,    11d17s23
     $     bc,ibc)                                                      11d17s23
      implicit real*8 (a-h,o-z)                                         11d17s23
      dimension idoubo(*),irefo(*),noc(*),nvirt(*),isblk1(4,*),i3x(*),  11d17s23
     $     itt(*),nbasdws(*),multh(8,8),vdinout(*),nfdat(5,4,*),        11d17s23
     $     itmp(4,2),itmp2(4,2)                                         11d17s23
      include "common.store"                                            11d17s23
      common/kmfind/invk1(2,8,8,8,2)
      data loopx/1000/
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      ibcoffo=ibcoff                                                    11d17s23
      loop=0
      nrootm=nroot-1                                                    11d17s23
      ivtmp=ibcoff                                                      11d16s23
      ibcoff=ivtmp+ndoub*nroot                                          11d17s23
      call enough('dorbd4x.vtmp',bc,ibc)                                11d17s23
      igoul=itt(3)+1
      igoulg=1910851
      igoulgg=igoulg
      ivoff=1                                                           11d16s23
      jvtmp=ivtmp                                                       11d16s23
      do isb=1,nsymb                                                    11d16s23
       isbv12=multh(isb,isymmrci)                                       11d16s23
       do isbv1=1,nsymb                                                 11d16s23
        isbv2=multh(isbv1,isbv12)                                       11d16s23
        if(isbv2.ge.isbv1)then                                          11d16s23
         if(isbv1.eq.isbv2)then                                         11d16s23
          nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         11d16s23
          nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                         11d16s23
          do i=0,nfdat(3,1,isb)-1                                       11d16s23
           do ir=0,nrootm                                               11d17s23
            do iv=0,nvirt(isbv1)-1                                      11d16s23
             itri=((iv*(iv+1))/2)+iv                                    11d17s23
             iad=jvtmp+ir+nroot*(i+nfdat(3,1,isb)*itri)                                        11d17
             bc(iad)=vdinout(ivoff+iv)*srh                                   11d17s23
            end do                                                      11d16s23
            ivoff=ivoff+nvirt(isbv1)                                    11d16s23
           end do                                                       11d16s23
          end do                                                        11d16s23
          isw=0                                                         11d17s23
         else                                                           11d16s23
          nvvt=nvirt(isbv1)*nvirt(isbv2)                                11d17s23
          nvvs=nvvt                                                     11d17s23
          isw=1                                                         11d17s23
         end if                                                         11d16s23
         ii=0
         nvv=nvvs                                                       11d17s23
         ips=+1                                                         11d17s23
         do l=1,4                                                       11d16s23
          do i=0,nfdat(3,l,isb)-1                                       11d16s23
           do ir=0,nrootm                                               11d17s23
            do iv2=0,nvirt(isbv2)-1                                     11d17s23
             iv1top=iv2+isw*(nvirt(isbv1)-iv2)-1                        11d17s23
             do iv1=0,iv1top                                            11d17s23
              irec=iv1+nvirt(isbv1)*iv2                                 11d17s23
              itrit=((iv2*(iv2+ips))/2)+iv1                              11d17s23
              itrif=((iv2*(iv2-1))/2)+iv1                               11d17s23
              itrit=itrit+isw*(irec-itrit)                              11d17s23
              itrif=itrif+isw*(irec-itrif)                              11d17s23
              iad=jvtmp+ir+nroot*(i+nfdat(3,l,isb)*itrit)               11d17s23
              bc(iad)=vdinout(ivoff+itrif)                                   11d17s23
             end do                                                     11d17s23
            end do                                                      11d17s23
            ivoff=ivoff+nvvt                                            11d17s23
           end do
          end do                                                        11d16s23
          jvtmp=jvtmp+nvv*nroot*nfdat(3,l,isb)                          11d17s23
          ips=-1                                                        11d17s23
          nvv=nvvt                                                      11d17s23
         end do                                                         11d16s23
        end if                                                          11d16s23
       end do                                                           11d16s23
      end do                                                            11d16s23
      mvoff=ivtmp                                                       11d17s23
      do isb=1,nsymb                                                    11d17s23
       isbv12=multh(isb,isymmrci)                                       11d17s23
       ivoff=mvoff                                                      11d17s23
       do isbv1=1,nsymb                                                 11d17s23
        isbv2=multh(isbv1,isbv12)                                       11d17s23
        if(isbv2.ge.isbv1)then                                          11d17s23
         if(isbv2.eq.isbv1)then                                         11d17s23
          nvvs=(nvirt(isbv1)*(nvirt(isbv1)+1))/2                         11d17s23
          nvvt=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         11d17s23
          isw=0                                                         11d17s23
         else                                                           11d17s23
          nvvt=nvirt(isbv1)*nvirt(isbv2)                                 11d17s23
          nvvs=nvvt
          isw=1                                                         11d17s23
         end if                                                         11d17s23
         ibc0=ibcoff                                                    11d17s23
         do l=1,4                                                       11d17s23
          itmp(l,1)=ibcoff                                              11d17s23
          nn=nvirt(isbv2)*nfdat(3,l,isb)*nroot*noc(isbv1)               11d17s23
          itmp2(l,1)=itmp(l,1)+nn                                       11d17s23
          itmp(l,2)=itmp2(l,1)+nn                                       11d17s23
          nn=nvirt(isbv1)*nfdat(3,l,isb)*nroot*noc(isbv2)               11d17s23
          itmp2(l,2)=itmp(l,2)+nn                                       11d17s23
          ibcoff=itmp2(l,2)+nn                                          11d17s23
         end do                                                         11d17s23
         call enough('dorbd4x.tmp2',bc,ibc)                             11d17s23
         do iz=ibc0,ibcoff-1                                            11d17s23
          bc(iz)=0d0                                                    11d17s23
         end do                                                         11d17s23
         jvoff=mvoff                                                    11d17s23
         do jsbv1=1,nsymb                                               11d17s23
          jsbv2=multh(jsbv1,isbv12)                                     11d17s23
          if(jsbv2.ge.jsbv1)then                                        11d17s23
           if(jsbv2.eq.jsbv1)then
            mvvs=(nvirt(jsbv1)*(nvirt(jsbv1)+1))/2                       11d17s23
            mvvt=(nvirt(jsbv1)*(nvirt(jsbv1)-1))/2                       11d17s23
           else                                                         11d17s23
            mvvs=nvirt(jsbv1)*nvirt(jsbv2)                               11d17s23
            mvvt=mvvs                                                   11d17s23
           end if                                                       11d17s23
           if(noc(isbv1).gt.0)then                                      11d21s23
            i2eu1=invk1(1,isbv2,jsbv2,isbv1,2)
            icase1=invk1(2,isbv2,jsbv2,isbv1,2)
            is1=i2eu1
            call ilimts(noc(isblk1(3,is1)),nvirt(isblk1(4,is1)),         11d17s23
     $           mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               11d17s23
            nhere=ih+1-il                                                11d17s23
           else                                                         11d21s23
            nhere=0                                                     11d21s23
           end if                                                       11d21s23
           if(nhere.gt.0)then                                           11d17s23
            if(isblk1(1,is1).eq.isblk1(2,is1))then                      11d17s23
             nrow3=(nvirt(isblk1(1,is1))*(nvirt(isblk1(1,is1))+1))/2    11d17s23
             isw3=0                                                     11d17s23
            else                                                        11d17s23
             nrow3=nvirt(isblk1(1,is1))*nvirt(isblk1(2,is1))            11d17s23
             isw3=1                                                     11d17s23
            end if                                                      11d17s23
            if(icase1.eq.1)then                                         11d20s23
             i10=i1s                                                    11d17s23
             i1n=noc(isblk1(3,is1))                                     11d17s23
             ii=i3x(is1)                                                11d17s23
             do i2=i2s,i2e                                              11d17s23
              i2m=i2-1                                                  11d17s23
              if(i2.eq.i2e)i1n=i1e                                      11d17s23
              do i1=i10,i1n                                             11d17s23
               i1m=i1-1                                                 11d17s23
               ips=+1
               mvv=mvvs                                                 11d17s23
               jjvoff=jvoff                                             11d17s23
               do l=1,4                                                 11d17s23
                iu2bot=(i2m-min(ips,0))*(1-isw)                         11d17s23
                nrn=nroot*nfdat(3,l,isb)                                11d17s23
                nrm=nrn-1                                               11d17s23
                do iu2=iu2bot,nvirt(jsbv2)-1                                 11d17s23
                 irec=i2m+nvirt(jsbv1)*iu2                              11d17s23
                 itri=((iu2*(iu2+ips))/2)+i2m                           11d17s23
                 icol=itri+isw*(irec-itri)                              11d17s23
                 iadv=jjvoff+nrn*icol                                   11d17s23
                 do iv2=0,nvirt(isbv2)-1                                11d17s23
                  irec=iv2+nvirt(isbv2)*iu2                             11d17s23
                  ix=max(iu2,iv2)                                       11d17s23
                  in=min(iu2,iv2)                                       11d17s23
                  itri=((ix*(ix+1))/2)+in                               11d17s23
                  itri0=itri
                  itri=itri+isw3*(irec-itri)+ii                         11d17s23
                  jtmp=itmp(l,1)+nrn*(iv2+nvirt(isbv2)*i1m)             11d17s23
                  do iii=0,nrm                                           11d17s23
                   bc(jtmp+iii)=bc(jtmp+iii)+bc(itri)*bc(iadv+iii)      11d20s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end do                                                  11d17s23
                ips=-1                                                  11d17s23
                jjvoff=jjvoff+nrn*mvv                                   11d17s23
                mvv=mvvt                                                11d17s23
               end do                                                   11d17s23
               ii=ii+nrow3                                              11d17s23
              end do                                                    11d17s23
              i10=1                                                     11d17s23
             end do                                                     11d17s23
            else if(icase1.eq.2)then                                    8d6s24
             i10=i1s                                                    11d17s23
             i1n=noc(isblk1(3,is1))                                     11d17s23
             ii=i3x(is1)                                                11d17s23
             do i2=i2s,i2e                                              11d17s23
              i2m=i2-1                                                  11d17s23
              if(i2.eq.i2e)i1n=i1e                                      11d17s23
              do i1=i10,i1n                                             11d17s23
               i1m=i1-1                                                 11d17s23
               ips=+1
               mvv=mvvs                                                 11d17s23
               jjvoff=jvoff                                             11d17s23
               do l=1,4                                                 11d17s23
                iu2bot=(i2m-min(ips,0))*(1-isw)                          11d17s23
                nrn=nroot*nfdat(3,l,isb)                                11d17s23
                nrm=nrn-1                                               11d17s23
                do iv2=0,nvirt(isbv2)-1                                 11d17s23
                 do iu2=iu2bot,nvirt(jsbv2)-1                                 11d17s23
                  irec=i2m+nvirt(jsbv1)*iu2                              11d17s23
                  itri=((iu2*(iu2+ips))/2)+i2m                           11d17s23
                  icol=itri+isw*(irec-itri)                              11d17s23
                  iadv=jjvoff+nrn*icol                                   11d17s23
                  irec=iu2+nvirt(jsbv2)*iv2+ii                          11d17s23
                  jtmp=itmp(l,1)+nrn*(iv2+nvirt(isbv2)*i1m)             11d17s23
                  do iii=0,nrm                                           11d17s23
                   bc(jtmp+iii)=bc(jtmp+iii)+bc(irec)*bc(iadv+iii)         11d17s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end do                                                  11d17s23
                ips=-1                                                  11d17s23
                jjvoff=jjvoff+nrn*mvv                                   11d17s23
                mvv=mvvt                                                11d17s23
               end do                                                   11d17s23
               ii=ii+nrow3                                              11d17s23
              end do                                                    11d17s23
              i10=1                                                     11d17s23
             end do                                                     11d17s23
            end if                                                      11d20s23
           end if                                                       11d20s23
           if(noc(isbv2).gt.0)then                                      11d21s23
            i2eu2=invk1(1,isbv1,jsbv2,isbv2,2)
            icase2=invk1(2,isbv1,jsbv2,isbv2,2)
            is1=i2eu2                                                    11d20s23
            call ilimts(noc(isblk1(3,is1)),nvirt(isblk1(4,is1)),         11d17s23
     $           mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               11d17s23
            nhere=ih+1-il                                                11d17s23
           else                                                         11d21s23
            nhere=0                                                     11d21s23
           end if                                                       11d21s23
           if(nhere.gt.0)then                                           11d17s23
            if(isblk1(1,is1).eq.isblk1(2,is1))then                      11d17s23
             nrow3=(nvirt(isblk1(1,is1))*(nvirt(isblk1(1,is1))+1))/2    11d17s23
             isw3=0                                                     11d17s23
            else                                                        11d17s23
             nrow3=nvirt(isblk1(1,is1))*nvirt(isblk1(2,is1))            11d17s23
             isw3=1                                                     11d17s23
            end if                                                      11d17s23
            if(icase2.eq.1)then                                         11d20s23
             i10=i1s                                                    11d17s23
             i1n=noc(isblk1(3,is1))                                     11d17s23
             ii=i3x(is1)                                                11d17s23
             do i2=i2s,i2e                                              11d17s23
              i2m=i2-1                                                  11d17s23
              if(i2.eq.i2e)i1n=i1e                                      11d17s23
              do i1=i10,i1n                                             11d17s23
               i1m=i1-1                                                 11d17s23
               ips=+1
               mvv=mvvs                                                 11d17s23
               jjvoff=jvoff                                             11d17s23
               do l=1,4                                                 11d17s23
                iu2bot=(i2m-min(ips,0))*(1-isw)                         11d17s23
                nrn=nroot*nfdat(3,l,isb)                                11d17s23
                nrm=nrn-1                                               11d17s23
                do iu2=iu2bot,nvirt(jsbv2)-1                                 11d17s23
                 irec=i2m+nvirt(jsbv1)*iu2                              11d17s23
                 itri=((iu2*(iu2+ips))/2)+i2m                           11d17s23
                 icol=itri+isw*(irec-itri)                              11d17s23
                 iadv=jjvoff+nrn*icol                                   11d17s23
                 do iv1=0,nvirt(isbv1)-1                                11d17s23
                  irec=iv1+nvirt(isbv1)*iu2                             11d17s23
                  ix=max(iu2,iv1)                                       11d17s23
                  in=min(iu2,iv1)                                       11d17s23
                  itri=((ix*(ix+1))/2)+in                               11d17s23
                  itri=itri+isw3*(irec-itri)+ii                         11d17s23
                  jtmp=itmp(l,2)+nrn*(iv1+nvirt(isbv1)*i1m)             11d17s23
                  do iii=0,nrm                                           11d17s23
                   bc(jtmp+iii)=bc(jtmp+iii)+bc(itri)*bc(iadv+iii)      11d20s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end do                                                  11d17s23
                ips=-1                                                  11d17s23
                jjvoff=jjvoff+nrn*mvv                                   11d17s23
                mvv=mvvt                                                11d17s23
               end do                                                   11d17s23
               ii=ii+nrow3                                              11d17s23
              end do                                                    11d17s23
              i10=1                                                     11d17s23
             end do                                                     11d17s23
            else if(icase2.eq.2)then                                    8d6s24
             i10=i1s                                                    11d17s23
             i1n=noc(isblk1(3,is1))                                     11d17s23
             ii=i3x(is1)                                                11d17s23
             do i2=i2s,i2e                                              11d17s23
              i2m=i2-1                                                  11d17s23
              if(i2.eq.i2e)i1n=i1e                                      11d17s23
              do i1=i10,i1n                                             11d17s23
               i1m=i1-1                                                 11d17s23
               ips=+1
               mvv=mvvs                                                 11d17s23
               jjvoff=jvoff                                             11d17s23
               do l=1,4                                                 11d17s23
                iu2bot=(i2m-min(ips,0))*(1-isw)                          11d17s23
                nrn=nroot*nfdat(3,l,isb)                                11d17s23
                nrm=nrn-1                                               11d17s23
                do iv1=0,nvirt(isbv1)-1                                 11d17s23
                 do iu2=iu2bot,nvirt(jsbv2)-1                                 11d17s23
                  irec=i2m+nvirt(jsbv1)*iu2                              11d17s23
                  itri=((iu2*(iu2+ips))/2)+i2m                           11d17s23
                  icol=itri+isw*(irec-itri)                              11d17s23
                  iadv=jjvoff+nrn*icol                                   11d17s23
                  irec=iu2+nvirt(jsbv2)*iv1+ii                          11d17s23
                  jtmp=itmp(l,2)+nrn*(iv1+nvirt(isbv1)*i1m)             11d17s23
                  do iii=0,nrm                                           11d17s23
                   bc(jtmp+iii)=bc(jtmp+iii)+bc(irec)*bc(iadv+iii)      11d20s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end do                                                  11d17s23
                ips=-1                                                  11d17s23
                jjvoff=jjvoff+nrn*mvv                                   11d17s23
                mvv=mvvt                                                11d17s23
               end do                                                   11d17s23
               ii=ii+nrow3                                              11d17s23
              end do                                                    11d17s23
              i10=1                                                     11d17s23
             end do                                                     11d17s23
            end if
           end if
           if(noc(isbv2).gt.0)then                                      11d21s23
            i2eu3=invk1(1,isbv1,jsbv1,isbv2,2)
            icase3=invk1(2,isbv1,jsbv1,isbv2,2)
            is1=i2eu3                                                    11d20s23
            call ilimts(noc(isblk1(3,is1)),nvirt(isblk1(4,is1)),         11d17s23
     $           mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               11d17s23
            nhere=ih+1-il                                                11d17s23
           else                                                         11d21s23
            nhere=0                                                     11d21s23
           end if                                                       11d21s23
           if(nhere.gt.0)then                                           11d17s23
            if(isblk1(1,is1).eq.isblk1(2,is1))then                      11d17s23
             nrow3=(nvirt(isblk1(1,is1))*(nvirt(isblk1(1,is1))+1))/2    11d17s23
             isw3=0                                                     11d17s23
            else                                                        11d17s23
             nrow3=nvirt(isblk1(1,is1))*nvirt(isblk1(2,is1))            11d17s23
             isw3=1                                                     11d17s23
            end if                                                      11d17s23
            if(icase3.eq.1)then                                         11d20s23
             i10=i1s                                                    11d17s23
             i1n=noc(isblk1(3,is1))                                     11d17s23
             ii=i3x(is1)                                                11d17s23
             do i2=i2s,i2e                                              11d17s23
              i2m=i2-1                                                  11d17s23
              if(i2.eq.i2e)i1n=i1e                                      11d17s23
              do i1=i10,i1n                                             11d17s23
               i1m=i1-1                                                 11d17s23
               ips=+1
               mvv=mvvs                                                 11d17s23
               jjvoff=jvoff                                             11d17s23
               do l=1,4                                                 11d17s23
                i2mm=i2m+min(ips,0)                                     11d21s23
                iu1top=i2mm+isw*(nvirt(jsbv1)-1-i2mm)                   11d21s23
                nrn=nroot*nfdat(3,l,isb)                                11d17s23
                nrm=nrn-1                                               11d17s23
                do iu1=0,iu1top                                         11d21s23
                 irec=iu1+nvirt(jsbv1)*i2m                              11d21s23
                 itri=((i2m*(i2m+ips))/2)+iu1                           11d21s23
                 icol=itri+isw*(irec-itri)                              11d17s23
                 iadv=jjvoff+nrn*icol                                   11d17s23
                 do iv1=0,nvirt(isbv1)-1                                11d17s23
                  irec=iv1+nvirt(isbv1)*iu1                             11d21s23
                  ix=max(iu1,iv1)                                       11d21s23
                  in=min(iu1,iv1)                                       11d21s23
                  itri=((ix*(ix+1))/2)+in                               11d17s23
                  itri=itri+isw3*(irec-itri)+ii                         11d17s23
                  jtmp=itmp2(l,2)+nrn*(iv1+nvirt(isbv1)*i1m)             11d17s23
                  do iii=0,nrm                                           11d17s23
                   bc(jtmp+iii)=bc(jtmp+iii)+bc(itri)*bc(iadv+iii)      11d20s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end do                                                  11d17s23
                ips=-1                                                  11d17s23
                jjvoff=jjvoff+nrn*mvv                                   11d17s23
                mvv=mvvt                                                11d17s23
               end do                                                   11d17s23
               ii=ii+nrow3                                              11d17s23
              end do                                                    11d17s23
              i10=1                                                     11d17s23
             end do                                                     11d17s23
            else if(icase3.eq.2)then                                    8d6s24
             i10=i1s                                                    11d17s23
             i1n=noc(isblk1(3,is1))                                     11d17s23
             ii=i3x(is1)                                                11d17s23
             do i2=i2s,i2e                                              11d17s23
              i2m=i2-1                                                  11d17s23
              if(i2.eq.i2e)i1n=i1e                                      11d17s23
              do i1=i10,i1n                                             11d17s23
               i1m=i1-1                                                 11d17s23
               ips=+1
               mvv=mvvs                                                 11d17s23
               jjvoff=jvoff                                             11d17s23
               do l=1,4                                                 11d17s23
                i2mm=i2m+min(ips,0)                                     11d21s23
                iu1top=i2mm+isw*(nvirt(jsbv1)-1-i2mm)                   11d21s23
                nrn=nroot*nfdat(3,l,isb)                                11d17s23
                nrm=nrn-1                                               11d17s23
                do iv1=0,nvirt(isbv1)-1                                 11d25s23
                 do iu1=0,iu1top                                         11d21s23
                  irec=iu1+nvirt(jsbv1)*i2m                              11d21s23
                  itri=((i2m*(i2m+ips))/2)+iu1                           11d21s23
                  icol=itri+isw*(irec-itri)                              11d17s23
                  iadv=jjvoff+nrn*icol                                   11d17s23
                  irec=iu1+nvirt(jsbv1)*iv1                             11d25s23
                  ix=max(iu1,iv1)                                       11d21s23
                  in=min(iu1,iv1)                                       11d21s23
                  itri=((ix*(ix+1))/2)+in                               11d17s23
                  itri=itri+isw3*(irec-itri)+ii                         11d17s23
                  jtmp=itmp2(l,2)+nrn*(iv1+nvirt(isbv1)*i1m)             11d17s23
                  do iii=0,nrm                                           11d17s23
                   bc(jtmp+iii)=bc(jtmp+iii)+bc(itri)*bc(iadv+iii)      11d20s23
                  end do                                                11d17s23
                 end do                                                 11d17s23
                end do                                                  11d17s23
                ips=-1                                                  11d17s23
                jjvoff=jjvoff+nrn*mvv                                   11d17s23
                mvv=mvvt                                                11d17s23
               end do                                                   11d17s23
               ii=ii+nrow3                                              11d17s23
              end do                                                    11d17s23
              i10=1                                                     11d17s23
             end do                                                     11d17s23
            end if
           end if                                                       11d20s23
           if(noc(isbv1).gt.0)then                                      11d21s23
            i2eu4=invk1(1,isbv2,jsbv1,isbv1,2)
            icase4=invk1(2,isbv2,jsbv1,isbv1,2)
            is1=i2eu4                                                    11d20s23
            call ilimts(noc(isblk1(3,is1)),nvirt(isblk1(4,is1)),         11d17s23
     $           mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)               11d17s23
            nhere=ih+1-il                                                11d17s23
            if(nhere.gt.0)then                                           11d17s23
             if(isblk1(1,is1).eq.isblk1(2,is1))then                      11d17s23
              nrow3=(nvirt(isblk1(1,is1))*(nvirt(isblk1(1,is1))+1))/2    11d17s23
              isw3=0                                                     11d17s23
             else                                                        11d17s23
              nrow3=nvirt(isblk1(1,is1))*nvirt(isblk1(2,is1))            11d17s23
              isw3=1                                                     11d17s23
             end if                                                      11d17s23
             if(icase4.eq.1)then                                         11d20s23
              i10=i1s                                                    11d17s23
              i1n=noc(isblk1(3,is1))                                     11d17s23
              ii=i3x(is1)                                                11d17s23
              do i2=i2s,i2e                                              11d17s23
               i2m=i2-1                                                  11d17s23
               if(i2.eq.i2e)i1n=i1e                                      11d17s23
               do i1=i10,i1n                                             11d17s23
                i1m=i1-1                                                 11d17s23
                ips=+1
                mvv=mvvs                                                 11d17s23
                jjvoff=jvoff                                             11d17s23
                do l=1,4                                                 11d17s23
                 i2mm=i2m+min(ips,0)                                     11d21s23
                 iu1top=i2mm+isw*(nvirt(jsbv1)-1-i2mm)                   11d21s23
                 nrn=nroot*nfdat(3,l,isb)                                11d17s23
                 nrm=nrn-1                                               11d17s23
                 do iu1=0,iu1top                                         11d21s23
                  irec=iu1+nvirt(jsbv1)*i2m                              11d21s23
                  itri=((i2m*(i2m+ips))/2)+iu1                           11d21s23
                  icol=itri+isw*(irec-itri)                              11d17s23
                  iadv=jjvoff+nrn*icol                                   11d17s23
                  do iv2=0,nvirt(isbv2)-1                                11d17s23
                   irec=iv2+nvirt(isbv2)*iu1                             11d21s23
                   ix=max(iu1,iv2)                                       11d21s23
                   in=min(iu1,iv2)                                       11d21s23
                   itri=((ix*(ix+1))/2)+in                               11d21s23
                   itri=itri+isw3*(irec-itri)+ii                         11d17s23
                   jtmp=itmp2(l,1)+nrn*(iv2+nvirt(isbv2)*i1m)            11d21s23
                   do iii=0,nrm                                          11d21s23
                    bc(jtmp+iii)=bc(jtmp+iii)+bc(itri)*bc(iadv+iii)      11d20s23
                   end do                                                11d17s23
                  end do                                                 11d17s23
                 end do                                                  11d17s23
                 ips=-1                                                  11d17s23
                 jjvoff=jjvoff+nrn*mvv                                   11d17s23
                 mvv=mvvt                                                11d17s23
                end do                                                   11d17s23
                ii=ii+nrow3                                              11d17s23
               end do                                                    11d17s23
               i10=1                                                     11d17s23
              end do                                                     11d17s23
             else if(icase4.eq.2)then                                   8d6s24
              i10=i1s                                                    11d17s23
              i1n=noc(isblk1(3,is1))                                     11d17s23
              ii=i3x(is1)                                                11d17s23
              do i2=i2s,i2e                                              11d17s23
               i2m=i2-1                                                  11d17s23
               if(i2.eq.i2e)i1n=i1e                                      11d17s23
               do i1=i10,i1n                                             11d17s23
                i1m=i1-1                                                 11d17s23
                ips=+1
                mvv=mvvs                                                 11d17s23
                jjvoff=jvoff                                             11d17s23
                do l=1,4                                                 11d17s23
                 i2mm=i2m+min(ips,0)                                     11d21s23
                 iu1top=i2mm+isw*(nvirt(jsbv1)-1-i2mm)                   11d21s23
                 nrn=nroot*nfdat(3,l,isb)                                11d17s23
                 nrm=nrn-1                                               11d17s23
                 do iv2=0,nvirt(isbv2)-1                                 11d21s23
                  do iu1=0,iu1top                                        11d21s23
                   irec=iu1+nvirt(jsbv1)*i2m                             11d21s23
                   itri=((i2m*(i2m+ips))/2)+iu1                          11d21s23
                   icol=itri+isw*(irec-itri)                              11d17s23
                   iadv=jjvoff+nrn*icol                                   11d17s23
                   irec=iu1+nvirt(jsbv1)*iv2+ii                          11d21s23
                   jtmp=itmp2(l,1)+nrn*(iv2+nvirt(isbv2)*i1m)            11d21s23
                   do iii=0,nrm                                           11d17s23
                    bc(jtmp+iii)=bc(jtmp+iii)+bc(irec)*bc(iadv+iii)         11d17s23
                   end do                                                11d17s23
                  end do                                                 11d17s23
                 end do                                                  11d17s23
                 ips=-1                                                  11d17s23
                 jjvoff=jjvoff+nrn*mvv                                   11d17s23
                 mvv=mvvt                                                11d17s23
                end do                                                   11d17s23
                ii=ii+nrow3                                              11d17s23
               end do                                                    11d17s23
               i10=1                                                     11d17s23
              end do                                                     11d17s23
             end if                                                      11d20s23
            end if                                                      11d21s23
           end if                                                       11d20s23
           mvv=mvvs                                                     11d17s23
           do l=1,4                                                     11d17s23
            nrn=nroot*nfdat(3,l,isb)                                    11d17s23
            jvoff=jvoff+nrn*mvv                                         11d17s23
            mvv=mvvt                                                    11d17s23
           end do                                                       11d17s23
          end if                                                        11d17s23
         end do                                                         11d17s23
         nvv=nvvs                                                       11d17s23
         ips=+1                                                         11d17s23
         do l=1,4                                                       11d17s23
          nrn=nroot*nfdat(3,l,isb)                                      11d17s23
          nrm=nrn-1                                                     11d17s23
          if(l.eq.1)then                                                11d17s23
           do i=0,nrn*nvirt(isbv2)*noc(isbv1)-1                         11d17s23
            bc(itmp(l,1)+i)=2d0*(bc(itmp(l,1)+i)+bc(itmp2(l,1)+i))      11d17s23
           end do                                                       11d17s23
           do i=0,nrn*nvirt(isbv1)*noc(isbv2)-1                         11d17s23
            bc(itmp(l,2)+i)=2d0*(bc(itmp(l,2)+i)+bc(itmp2(l,2)+i))      11d17s23
           end do                                                       11d17s23
          else                                                          11d17s23
           do i=0,nrn*nvirt(isbv2)*noc(isbv1)-1                         11d17s23
            bc(itmp(l,1)+i)=2d0*(bc(itmp(l,1)+i)-bc(itmp2(l,1)+i))      11d17s23
           end do                                                       11d17s23
           do i=0,nrn*nvirt(isbv1)*noc(isbv2)-1                         11d17s23
            bc(itmp(l,2)+i)=2d0*(bc(itmp2(l,2)+i)-bc(itmp(l,2)+i))      11d17s23
           end do                                                       11d17s23
          end if                                                        11d17s23
          do iv2=0,nvirt(isbv2)-1                                       11d17s23
           iv2p=iv2+noc(isbv2)                                          11d20s23
           iv2x=iv2+min(ips,0)                                          11d17s23
           iv1top=iv2x+isw*(nvirt(isbv1)-1-iv2x)                        11d17s23
           do iv1=0,iv1top                                              11d17s23
            iv1p=iv1+noc(isbv1)                                         11d20s23
            irec=iv1+nvirt(isbv1)*iv2                                   11d17s23
            itri=((iv2*(iv2+ips))/2)+iv1                                11d17s23
            itri=itri+isw*(irec-itri)                                   11d17s23
            iadv=ivoff+nrn*itri                                         11d17s23
            do ie=0,noc(isbv1)-1                                        11d17s23
             iad1=itt(isbv1)+iv1p+nbasdws(isbv1)*ie                     11d20s23
             iadt=itmp(l,1)+nrn*(iv2+nvirt(isbv2)*ie)                   11d17s23
             do i=0,nrm                                                 11d17s23
              bc(iad1)=bc(iad1)+bc(iadv+i)*bc(iadt+i)                   11d17s23
             end do                                                     11d17s23
            end do                                                      11d17s23
            do ie=0,noc(isbv2)-1                                        11d17s23
             iad1=itt(isbv2)+iv2p+nbasdws(isbv2)*ie                     11d20s23
             iadt=itmp(l,2)+nrn*(iv1+nvirt(isbv1)*ie)                   11d17s23
             do i=0,nrm                                                 11d17s23
              bc(iad1)=bc(iad1)+bc(iadv+i)*bc(iadt+i)                   11d17s23
             end do                                                     11d17s23
            end do                                                      11d17s23
           end do                                                       11d17s23
          end do                                                        11d17s23
          ivoff=ivoff+nrn*nvv                                           11d17s23
          nvv=nvvt                                                      11d17s23
          ips=-1                                                        11d17s23
         end do                                                         11d17s23
         ibcoff=ibc0                                                    11d17s23
         mvoffn=ivoff                                                   11d21s23
        end if                                                          11d17s23
       end do                                                           11d17s23
       mvoff=mvoffn                                                     11d21s23
      end do                                                            11d17s23
      ibcoff=ibcoffo                                                    11d17s23
      return                                                            11d17s23
      end                                                               11d17s23
