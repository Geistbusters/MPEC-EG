c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hddcsf12b(nhand,ih0av,nh0av,iff2,nff2,mdon,            6d13s21
     $     mdoo,ioooo,jmats,kmats,nec,nvirt,shift,nsymb,multh,bc,ibc)   11d10s22
      implicit real*8 (a-h,o-z)                                         7d15s19
      external second                                                   11d1s19
      logical ldebug                                                    5d28s21
      integer*8 nhand(mdoo+1,*),i1,i12,in,i1off                                6d8s21
      dimension ih0av(*),nh0av(*),ioooo(*),jmats(*),kmats(*),           5d28s21
     $     iff2(*),nff2(mdoo+1,*),idorb(32),isorb(32),nvirt(*),         6d13s21
     $     multh(8,8)                                                   6d13s21
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      common/kmfind/invk1(2,8,8,8,2)                                    6d30s18
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      include "common.mrci"                                             7d15s19
      include "common.store"                                            7d15s19
      i1=1                                                              11d19s20
      iffoff=1                                                          6d13s21
      do isb=1,nsymb                                                    6d13s21
       isbv12=multh(isb,isymmrci)                                       6d13s21
       nvisv=0                                                          6d13s21
       nvnotv=0                                                         6d13s21
       do isbv1=1,nsymb                                                 6d13s21
        isbv2=multh(isbv1,isbv12)                                       6d13s21
        if(isbv2.ge.isbv1)then                                          6d13s21
         if(isbv1.eq.isbv2)then                                         6d13s21
          nvisv=nvisv+nvirt(isbv1)                                      6d13s21
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         6d13s21
         else                                                           6d13s21
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 6d13s21
         end if                                                         6d13s21
         nvnotv=nvnotv+nvv                                              6d13s21
        end if                                                          6d13s21
       end do                                                           6d13s21
       nall=nvisv+nvnotv                                                6d13s21
       do nclop=1,mdoo+1                                                 11d19s20
        nclo=nclop-1                                                     11d19s20
        if(nff2(nclop,isb).gt.0)then                                    6d13s21
         igg=ibcoff                                                        6d13s21
         ibcoff=igg+nall*nff2(nclop,isb)                                6d13s21
         call enough('hddcsf12b.  1',bc,ibc)
         nopen=nec-nclo*2-2                                              11d19s20
         iarg=nclop-mdon                                                 11d19s20
         do i=igg,ibcoff-1                                              6d13s21
          bc(i)=0d0                                                     6d13s21
         end do                                                         6d13s21
         jgg=igg                                                         11d19s20
         do if=1,nff2(nclop,isb)                                        6d13s21
          ncheck=popcnt(iff2(iffoff))                                    5d28s21
          if(ncheck.ne.nclo)stop 'ncheckc'
          ii=1                                                           11d19s20
          do i=1,norb                                                    11d19s20
           if(btest(iff2(iffoff),i))then                                 5d28s21
            idorb(ii)=i                                                  11d19s20
            ii=ii+1                                                      11d19s20
           end if                                                        11d19s20
          end do                                                         11d19s20
          iffoff=iffoff+1                                                11d19s20
          ii=1                                                           11d19s20
          ncheck=popcnt(iff2(iffoff))                                    5d28s21
          if(ncheck.ne.nopen)stop 'nhecko'
           do i=1,norb                                                    11d19s20
            if(btest(iff2(iffoff),i))then                                 5d28s21
            isorb(ii)=i                                                  11d19s20
            ii=ii+1                                                      11d19s20
           end if                                                        11d19s20
          end do                                                         11d19s20
          iffoff=iffoff+1                                                11d19s20
          sumc=0d0                                                         7d15s19
          do i=1,nclo                                                    11d19s20
           ig=irel(idorb(i))-1                                           11d19s20
           is=ism(idorb(i))                                              11d19s20
           iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
           h0=bc(iad)                                                      7d15s19
           do j=1,nclo                                                   11d19s20
            jg=irel(idorb(j))-1                                          11d19s20
            js=ism(idorb(j))                                             11d19s20
            xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)             11d15s22
            xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)             11d15s22
            sumc=sumc+2d0*xj-xk                                            7d15s19
           end do                                                          7d15s19
           sumc=sumc+h0*2d0                                                7d15s19
          end do                                                           7d15s19
          sumo=0d0                                                         7d15s19
          do i=1,nopen                                                   11d19s20
           ig=irel(isorb(i))-1                                           11d19s20
           is=ism(isorb(i))                                              11d19s20
           iad=ih0av(is)+ig*(nh0av(is)+1)                                  7d15s19
           h0=bc(iad)                                                      7d15s19
           sumo=sumo+h0                                                    7d15s19
           do j=1,nopen                                                  11d19s20
            if(i.ne.j)then                                                 7d15s19
             jg=irel(isorb(j))-1                                         11d19s20
             js=ism(isorb(j))                                            11d19s20
             xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)            11d15s22
             sumo=sumo+xj*0.5d0                                            7d15s19
            end if                                                         7d15s19
           end do                                                          7d15s19
          end do                                                           7d15s19
          sum=sumc+sumo+shift                                              7d15s19
          do i=1,nclo                                                    11d19s20
           is=ism(idorb(i))                                              11d19s20
           ig=irel(idorb(i))-1                                           11d19s20
           do j=1,nopen                                                  11d19s20
            js=ism(isorb(j))                                             11d19s20
            jg=irel(isorb(j))-1                                          11d19s20
            xj=getint(ioooo,is,is,js,js,ig,ig,jg,jg,bc,ibc)             11d15s22
            xk=getint(ioooo,is,js,js,is,ig,jg,jg,ig,bc,ibc)             11d15s22
            sum=sum+xj*2d0-xk                                              7d15s19
           end do                                                          7d15s19
          end do                                                           7d15s19
          do isbv1=1,nsymb                                              6d13s21
           isbv2=multh(isbv1,isbv12)                                    6d13s21
           if(isbv2.ge.isbv1)then                                       6d13s21
            call ilimts(nvirt(isbv1),nvirt(isbv1),mynprocg,mynowprog,   6d13s21
     $          il,ih,i1s,i1e,i2s,i2e)                                  6d13s21
            nhere=ih+1-il                                                     7d15s19
            call ilimts(nvirt(isbv2),nvirt(isbv2),mynprocg,mynowprog,
     $           jl,jh,i1s,i1e,i2s,i2e)                                 6d13s21
            nherej=jh+1-jl                                                     7d15s19
            if(isbv1.eq.isbv2)then                                          5d28s21
             nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                      6d13s21
             isw=0                                                            5d28s21
             do iv=0,nvirt(isbv1)-1                                     6d13s21
              icol=iv*(nvirt(isbv1)+1)+1                                6d13s21
              if(icol.ge.il.and.icol.le.ih)then                            5d28s21
               ivp=iv+irefo(isbv1)                                         5d28s21
               icolp=icol-il                                                  7d15s19
               iad=ih0av(isbv1)+ivp*(nh0av(isbv1)+1)                       5d28s21
               h0=bc(iad)                                                     7d15s19
               here=h0*2d0+sum                                             5d28s21
               do j=1,nclo                                                  11d19s20
                jg=irel(idorb(j))-1                                         11d19s20
                js=ism(idorb(j))                                            11d19s20
                i2eu=inv(1,js,js,isbv1)                                        7d15s19
                irow=((jg*(jg+1))/2)+jg                                       7d15s19
                iadj=jmats(i2eu)+icolp+nhere*irow                             7d15s19
                i2eu=invk1(1,js,js,isbv1,1)                                    7d15s19
                irow=jg+irefo(js)*jg                                          7d15s19
                iadk=kmats(i2eu)+icolp+nhere*irow                             7d15s19
                here=here+2d0*(2d0*bc(iadj)-bc(iadk))                      6d7s21
               end do                                                         7d15s19
               do j=1,nopen                                                5d28s21
                jg=irel(isorb(j))-1                                         11d19s20
                js=ism(isorb(j))                                            11d19s20
                i2eu=inv(1,js,js,isbv1)                                        7d15s19
                irow=((jg*(jg+1))/2)+jg                                       7d15s19
                iadj=jmats(i2eu)+icolp+nhere*irow                             7d15s19
                i2eu=invk1(1,js,js,isbv1,1)                                    7d15s19
                irow=jg+irefo(js)*jg                                          7d15s19
                iadk=kmats(i2eu)+icolp+nhere*irow                             7d15s19
                here=here+2d0*bc(iadj)-bc(iadk)                               7d15s19
               end do                                                      5d28s21
               bc(jgg+iv)=here                                             5d28s21
              end if                                                       5d28s21
             end do                                                        5d28s21
             jgg=jgg+nvirt(isbv1)                                       6d13s21
            else                                                            5d28s21
             nvv=nvirt(isbv1)*nvirt(isbv2)                              6d13s21
             isw=1                                                            5d28s21
            end if                                                          5d28s21
            do iv2=0,nvirt(isbv2)-1                                     6d13s21
             if(isbv1.eq.isbv2)then                                        5d28s21
              itop=iv2-1                                                   5d28s21
             else                                                          5d28s21
              itop=nvirt(isbv1)-1                                       6d13s21
             end if                                                        5d28s21
             iv2p=iv2+irefo(isbv2)                                         5d28s21
             iad2=ih0av(isbv2)+iv2p*(nh0av(isbv2)+1)                       5d28s21
             do iv1=0,itop                                                 5d28s21
              irec=iv1+nvirt(isbv1)*iv2                                          5d28s21
              itri=((iv2*(iv2-1))/2)+iv1                                   5d28s21
              lgg=jgg+itri+isw*(irec-itri)                                 5d28s21
              icol=iv1+nvirt(isbv1)*iv1+1                                        5d28s21
              here=0d0                                                     5d28s21
              if(icol.ge.il.and.icol.le.ih)then                            5d28s21
               icolp=icol-il                                               5d28s21
               iv1p=iv1+irefo(isbv1)                                       5d28s21
               iad1=ih0av(isbv1)+iv1p*(nh0av(isbv1)+1)                     5d28s21
               here=bc(iad1)+bc(iad2)+sum                                  5d28s21
               do j=1,nopen                                                5d28s21
                jg=irel(isorb(j))-1                                         11d19s20
                js=ism(isorb(j))                                            11d19s20
                i2eu=inv(1,js,js,isbv1)                                         7d15s19
                irow=((jg*(jg+1))/2)+jg                                       7d15s19
                iad=jmats(i2eu)+icolp+nhere*irow                              7d15s19
                here=here+bc(iad)                                             7d15s19
               end do                                                      5d28s21
               do j=1,nclo                                                  11d19s20
                jg=irel(idorb(j))-1                                         11d19s20
                js=ism(idorb(j))                                            11d19s20
                i2eu=inv(1,js,js,isbv1)                                        7d15s19
                irow=((jg*(jg+1))/2)+jg                                       7d15s19
                iadj=jmats(i2eu)+icolp+nhere*irow                             7d15s19
                i2eu=invk1(1,js,js,isbv1,1)                                    7d15s19
                irow=jg+irefo(js)*jg                                          7d15s19
                iadk=kmats(i2eu)+icolp+nhere*irow                             7d15s19
                here=here+2d0*bc(iadj)-bc(iadk)                               7d15s19
               end do                                                         7d15s19
              end if                                                       5d28s21
              icol=iv2+nvirt(isbv2)*iv2+1                               6d13s21
              if(icol.ge.jl.and.icol.le.jh)then                            5d28s21
               icolp=icol-jl                                               5d28s21
               do j=1,nopen                                                5d28s21
                jg=irel(isorb(j))-1                                         11d19s20
                js=ism(isorb(j))                                            11d19s20
                i2eu=inv(1,js,js,isbv2)                                         7d15s19
                irow=((jg*(jg+1))/2)+jg                                       7d15s19
                iad=jmats(i2eu)+icolp+nherej*irow                              7d15s19
                here=here+bc(iad)                                             7d15s19
               end do                                                      5d28s21
               do j=1,nclo                                                  11d19s20
                jg=irel(idorb(j))-1                                         11d19s20
                js=ism(idorb(j))                                            11d19s20
                i2eu=inv(1,js,js,isbv2)                                        7d15s19
                irow=((jg*(jg+1))/2)+jg                                       7d15s19
                iadj=jmats(i2eu)+icolp+nherej*irow                             7d15s19
                i2eu=invk1(1,js,js,isbv2,1)                                    7d15s19
                irow=jg+irefo(js)*jg                                          7d15s19
                iadk=kmats(i2eu)+icolp+nherej*irow                             7d15s19
                here=here+2d0*bc(iadj)-bc(iadk)                               7d15s19
               end do                                                         7d15s19
              end if                                                       5d28s21
              bc(lgg)=here                                                 5d28s21
             end do                                                        5d28s21
            end do                                                         5d28s21
            jgg=jgg+nvv                                                 6d13s21
           end if                                                       6d13s21
          end do                                                        6d13s21
         end do                                                          11d19s20
         in=nff2(nclop,isb)                                             6d13s21
         i12=nall                                                       6d13s21
         call ddi_acc(bc,ibc,nhand(nclop,isb),i1,i12,i1,in,bc(igg))     11d15s22
         ibcoff=igg                                                     6d13s21
        end if                                                           11d19s20
       end do                                                            11d19s20
      end do                                                            6d13s21
      return                                                            7d15s19
      end                                                               7d15s19
