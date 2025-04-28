c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine tt4v(itt4v,iaddr,naddr,eraw2,nsymb,multh,nbasisp,      5d20s24
     $           nerimat,lsa,lsb,lsc,lsd,nan,iptap,nbn,iptbp,in,        5d31s24
     $           ncomp,lsame,noc,nvirt,bc,ibc)                          6d5s24
      implicit real*8 (a-h,o-z)                                         5d20s24
      external second                                                   8d23s24
      include "common.store"
      logical lsame,ldebug,look                                              5d31s24
      dimension itt4v(8,*),iaddr(36,2,*),naddr(36,*),eraw2(*),          8d23s24
     $     multh(8,8),nbasisp(*),iptap(*),iptbp(*),nan(*),nbn(*),       6d7s24
     $     nvirt(*),noc(*)                                              6d7s24
      common/timerocm/tovrx,telapo(15)                                  8d23s24
c
c     integrals are nan,nbn,nbasisp(lsc),nbasdws(lsd)
      data loop,loopx/0,36/
      data icall/0/
      save
      icall=icall+1
      look=.false.
      igoal=1987803
      igoal2=2212059
      ioff=0
      ifound=0
c
c     assuming order is LL,LS,SL,SS in vectors
      if(in.eq.0)then                                                   5d21s24
       ier2=1                                                           5d29s24
       iqshift=0                                                        5d29s24
      else if(in.eq.1)then                                              5d21s24
       ier2=2                                                           5d29s24
       iqshift=0                                                        5d29s24
      else if(in.eq.2)then                                              5d21s24
       ier2=3                                                           5d29s24
       iqshift=1                                                        5d29s24
      else                                                              5d21s24
       ier2=4                                                           5d29s24
       iqshift=1                                                        5d29s24
      end if                                                            5d21s24
      look=ifound.ne.0
cc    j1
      isbv1=lsa                                                         5d31s24
      isbv2=lsd                                                         5d31s24
      jsbv1=lsb                                                         5d31s24
      jsbv2=lsc                                                         5d31s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d6s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
       if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
        nv1=nvirt(isbv1)                                                6d5s24
        nv2=nvirt(isbv2)                                                6d5s24
        mv1=nbasisp(jsbv1)                                              6d10s24
        mv2=nbasisp(jsbv2)                                              6d10s24
        if(isbv2.eq.isbv1)then                                          6d7s24
         isw=0                                                          5d20s24
         nvvs=(nv1*(nv1+1))/2                                           5d31s24
         nvvt=(nv1*(nv1-1))/2                                           5d31s24
        else                                                            5d20s24
         isw=1                                                          5d20s24
         nvvs=nv1*nv2                                                   5d31s24
         nvvt=nvvs                                                      5d20s24
        end if                                                          5d20s24
        isub=1                                                          5d31s24
        nvv=nvvs                                                        5d31s24
        mvv=mv1*mv2                                                     6d7s24
        do ist=1,2                                                      5d31s24
         if(min(nvv,naddr(jtv12,ist)).gt.0)then                         8d26s24
          if(itt4v(lsa,1).gt.0)then                                     8d26s24
          call second(time1)
          do jj=0,nbn(lsb)-1                                            5d31s24
           jv1=ibc(iptbp(lsb)+jj)                                       5d31s24
           do jv2=0,nbasisp(jsbv2)-1                                    6d7s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do iv2=0,nvirt(isbv2)-1                                     6d6s24
             iv1top=iv2+1-ist                                           6d5s24
             iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
             iv2p=iv2+noc(isbv2)                                        6d5s24
             do iv1=0,iv1top                                            6d5s24
              irec=iv1+nvirt(isbv1)*iv2                                 6d5s24
              itri=((iv2*(iv2+isub))/2)+iv1                             6d5s24
              itri=itri+isw*(irec-itri)                                 6d5s24
              iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri              6d11s24
              sum=0d0                                                   6d11s24
              do i=0,naddr(jtv12,ist)-1                                 6d11s24
               orig=sum
               sum=sum+bc(jad+i)*bc(iad+i)                              6d11s24
              end do                                                    6d11s24
              do ii=0,nan(lsa)-1                                        5d31s24
               iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)               5d31s24
               iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv2                     5d31s24
     $                    +nbasisp(jsbv2)*iv2p))                        6d5s24
               iadtt=itt4v(lsa,1)+iv1+nvirt(isbv1)*iq                   6d7s24
               bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                      5d31s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          call second(time2)
          telapo(1)=telapo(1)+time2-time1-tovr
          else                                                          8d27s24
          call second(time2)
          end if                                                        8d26s24
          if(itt4v(lsa,5).gt.0)then                                     8d27s24
           nrow=nbn(lsb)*nbasisp(jsbv2)                                 8d26s24
           itmp=ibcoff                                                  8d26s24
           isum=itmp+nrow*naddr(jtv12,ist)                              8d26s24
           ibcoff=isum+nrow*nvv                                         8d26s24
           call enough('tt4v.sum12b',bc,ibc)                            8d26s24
           jtmp=itmp                                                    8d26s24
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv1=ibc(iptbp(lsb)+jj)                                       5d31s24
            do jv2=0,nbasisp(jsbv2)-1                                    6d7s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
             do i=0,naddr(jtv12,ist)-1                                  8d26s24
              bc(jtmp+i)=bc(jad+i)                                      8d26s24
             end do                                                     8d26s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d26s24
            end do                                                      8d26s24
           end do                                                       8d26s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),1d0,            8d26s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d26s24
     $         bc(isum),nvv,'tt4v.sum12b')                              8d26s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv2)                       8d26s24
           do ii=0,nan(lsa)-1                                           8d23s24
            iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)                  8d23s24
            iadtt=itt4v(lsa,5)+nvirt(isbv1)*iq                          8d23s24
            do jj=0,nbn(lsb)-1                                            5d31s24
             do jv2=0,nbasisp(jsbv2)-1                                    6d7s24
              lsum=isum+nvv*(jv2+nbasisp(jsbv2)*jj)                     8d26s24
              iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv2                      8d23s24
     $                    +nbasisp(jsbv2)*noc(isbv2)))                        6d5s24
              do iv2=0,nvirt(isbv2)-1                                     6d6s24
               iv1top=iv2+1-ist                                           6d5s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                8d27s24
               jadv=iadv+iv2*nrowi
               irec=nvirt(isbv1)*iv2                                    8d27s24
               itri=((iv2*(iv2+isub))/2)                                8d27s24
               itri=itri+isw*(irec-itri)+lsum                           8d27s24
               do iv1=0,iv1top                                            6d5s24
                bc(iadtt+iv1)=bc(iadtt+iv1)+bc(itri+iv1)*eraw2(jadv)     8d23s24
               end do                                                   8d26s24
              end do                                                    8d26s24
             end do                                                     8d26s24
            end do                                                      8d26s24
           end do                                                       8d26s24
           ibcoff=itmp                                                  8d26s24
           call second(time3)
           telapo(5)=telapo(5)+time3-time2-tovr
          end if                                                        8d26s24
         end if                                                         5d31s24
         isub=-1                                                        5d31s24
         nvv=nvvt                                                       5d31s24
        end do                                                          5d31s24
       end if                                                           5d31s24
      end if                                                            5d31s24
      if(.not.lsame)then                                                5d31s24
       isbv1=lsb                                                        5d31s24
       isbv2=lsd                                                        5d31s24
       jsbv1=lsa                                                        5d31s24
       jsbv2=lsc                                                        5d31s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         5d31s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
        if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
         nv1=nvirt(isbv1)                                               6d5s24
         nv2=nvirt(isbv2)                                               6d5s24
         mv1=nbasisp(jsbv1)                                             6d10s24
         mv2=nbasisp(jsbv2)                                             6d10s24
         if(isbv2.eq.isbv1)then                                         6d7s24
          isw=0                                                          5d20s24
          nvvs=(nv1*(nv1+1))/2                                          5d31s24
          nvvt=(nv1*(nv1-1))/2                                          5d31s24
         else                                                            5d20s24
          isw=1                                                          5d20s24
          nvvs=nv1*nv2                                                  5d31s24
          nvvt=nvvs                                                      5d20s24
         end if                                                          5d20s24
         isub=1                                                          5d31s24
         nvv=nvvs                                                        5d31s24
         mvv=mv1*mv2                                                    6d7s24
         nrow=nan(lsa)*nbasisp(jsbv2)                                   8d26s24
         do ist=1,2                                                      5d31s24
          if(min(nvv,naddr(jtv12,ist)).gt.0)then                        8d26s24
           call second(time1)
           if(itt4v(lsb,1).gt.0)then                                    8d26s24
           do jj=0,nan(lsa)-1                                            5d31s24
            jv1=ibc(iptap(lsa)+jj)                                       5d31s24
            do jv2=0,nbasisp(jsbv2)-1                                   6d7s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do iv2=0,nvirt(isbv2)-1                                    6d5s24
              iv2p=iv2+noc(isbv2)                                       6d5s24
              iv1top=iv2+1-ist                                           5d31s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
              do iv1=0,iv1top                                            5d31s24
               irec=iv1+nvirt(isbv1)*iv2                                6d5s24
               itri=((iv2*(iv2+isub))/2)+iv1                             5d31s24
               itri=itri+isw*(irec-itri)                                 5d31s24
               sum=0d0                                                  5d31s24
               iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri             6d11s24
               do i=0,naddr(jtv12,ist)-1                                5d31s24
                sum=sum+bc(jad+i)*bc(iad+i)                             6d11s24
               end do                                                   5d31s24
               do ii=0,nbn(lsb)-1                                       5d31s24
                iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)              5d31s24
                iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv2                    5d31s24
     $                    +nbasisp(jsbv2)*iv2p))                           5d21s24
                iadtt=itt4v(lsb,1)+iv1+nvirt(isbv1)*iq                  6d7s24
                bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                     5d31s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           call second(time2)
           telapo(1)=telapo(1)+time2-time1-tovr
           else                                                         8d26s24
            call second(time2)
           end if                                                       8d26s24
           if(itt4v(lsb,5).gt.0)then                                    8d27s24
           itmp=ibcoff                                                  8d26s24
           isum=itmp+nrow*naddr(jtv12,ist)                              8d26s24
           ibcoff=isum+nrow*nvv                                         8d26s24
           call enough('tt4v.Sum12',bc,ibc)                             8d26s24
           jtmp=itmp                                                    8d26s24
           do jj=0,nan(lsa)-1                                            5d31s24
            jv1=ibc(iptap(lsa)+jj)                                       5d31s24
            do jv2=0,nbasisp(jsbv2)-1                                   6d7s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do i=0,naddr(jtv12,ist)-1                                  8d26s24
              bc(jtmp+i)=bc(jad+i)                                      8d26s24
             end do                                                     8d26s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d26s24
            end do                                                      8d26s24
           end do                                                       8d26s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),1d0,            8d26s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d26s24
     $         bc(isum),nvv,'tt4v.Sum12')                               8d26s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv2)                       8d27s24
           do ii=0,nbn(lsb)-1                                           8d26s24
            iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)                  8d26s24
            iadtt=itt4v(lsb,5)+nvirt(isbv1)*iq                          8d26s24
            do jj=0,nan(lsa)-1                                            5d31s24
             jv1=ibc(iptap(lsa)+jj)                                       5d31s24
             do jv2=0,nbasisp(jsbv2)-1                                   6d7s24
              lsum=isum+nvv*(jv2+nbasisp(jsbv2)*jj)                     8d26s24
              iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv2                      8d27s24
     $                    +nbasisp(jsbv2)*noc(isbv2)))                  8d27s24
              do iv2=0,nvirt(isbv2)-1                                    6d5s24
               jadv=iadv+iv2*nrowi                                      8d27s24
               iv1top=iv2+1-ist                                           5d31s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
               irec=nvirt(isbv1)*iv2                                    8d26s24
               itri=((iv2*(iv2+isub))/2)                                8d26s24
               itri=itri+isw*(irec-itri)+lsum                           8d26s24
               do iv1=0,iv1top                                          8d23s24
                bc(iadtt+iv1)=bc(iadtt+iv1)+bc(itri+iv1)*eraw2(jadv)    8d27s24
               end do                                                   8d23s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           ibcoff=itmp                                                  8d26s24
           end if
           call second(time3)
           telapo(5)=telapo(5)+time3-time2-tovr
          end if                                                         5d31s24
          isub=-1                                                       5d31s24
          nvv=nvvt                                                      5d31s24
         end do                                                          5d31s24
        end if                                                           5d31s24
       end if                                                            5d31s24
      end if                                                            5d31s24
cc    j2
      isbv1=lsb                                                         5d31s24
      isbv2=lsd                                                         5d31s24
      jsbv1=lsa                                                         5d31s24
      jsbv2=lsc                                                         5d31s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d6s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
       if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
        nv1=nvirt(isbv1)                                                6d5s24
        nv2=nvirt(isbv2)                                                6d5s24
        mv1=nbasisp(jsbv1)                                              6d10s24
        mv2=nbasisp(jsbv2)                                              6d10s24
        if(isbv2.eq.isbv1)then                                          6d7s24
         isw=0                                                          5d20s24
         nvvs=(nv1*(nv1+1))/2                                           5d31s24
         nvvt=(nv1*(nv1-1))/2                                           5d31s24
        else                                                            5d20s24
         isw=1                                                          5d20s24
         nvvs=nv1*nv2                                                   5d31s24
         nvvt=nvvs                                                      5d20s24
        end if                                                          5d20s24
        isub=1                                                          5d31s24
        nvv=nvvs                                                        5d31s24
        mvv=mv1*mv2                                                     6d7s24
        nrow=nan(lsa)*nbasisp(jsbv2)                                    8d26s24
        ntogether=naddr(jtv12,1)+naddr(jtv12,2)                         8d26s24
        do ist=1,2                                                      5d31s24
         if(min(nvv,naddr(jtv12,ist)).gt.0)then                         8d26s24
          if(itt4v(lsb,2).gt.0)then                                     8d26s24
          call second(time1)
          do jj=0,nan(lsa)-1                                            5d31s24
           jv1=ibc(iptap(lsa)+jj)                                       5d31s24
           do jv2=0,nbasisp(jsbv2)-1                                    6d7s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do iv2=0,nvirt(isbv2)-1                                     6d6s24
             iv1top=iv2+1-ist                                           6d5s24
             iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
             iv2p=iv2+noc(isbv2)                                        6d5s24
             do iv1=0,iv1top                                            6d5s24
              irec=iv1+nvirt(isbv1)*iv2                                 6d5s24
              itri=((iv2*(iv2+isub))/2)+iv1                             6d5s24
              itri=itri+isw*(irec-itri)                                 6d5s24
              sum=0d0                                                   5d31s24
              iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri              6d11s24
              do i=0,naddr(jtv12,ist)-1                                 5d31s24
               sum=sum+bc(jad+i)*bc(iad+i)                              6d11s24
              end do                                                    5d31s24
              do ii=0,nbn(lsb)-1                                        5d31s24
               iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)               5d31s24
               iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv2                     5d31s24
     $                    +nbasisp(jsbv2)*iv2p))                        6d5s24
               iadtt=itt4v(lsb,2)+iv1+nvirt(isbv1)*iq                   6d7s24
               bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                      5d31s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          call second(time2)
          telapo(2)=telapo(2)+time2-time1-tovr
          else                                                          8d27s24
           call second(time2)
          end if                                                        8d26s24
          if(itt4v(lsb,6).gt.0)then                                     8d27s24
          itmp=ibcoff                                                   8d26s24
          isum=itmp+nrow*naddr(jtv12,ist)                               8d26s24
          ibcoff=isum+nvv*nrow                                          8d26s24
          call enough('tt4v.sum21',bc,ibc)                              8d23s24
          jtmp=itmp                                                     8d26s24
          do jj=0,nan(lsa)-1                                            5d31s24
           jv1=ibc(iptap(lsa)+jj)                                       5d31s24
           do jv2=0,nbasisp(jsbv2)-1                                    6d7s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do i=0,naddr(jtv12,ist)-1                                   8d26s24
             bc(jtmp+i)=bc(jad+i)                                       8d26s24
            end do                                                      8d26s24
            jtmp=jtmp+naddr(jtv12,ist)                                  8d26s24
           end do                                                       8d26s24
          end do                                                        8d26s24
          call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),1d0,             8d26s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d26s24
     $         bc(isum),nvv,'tt4v.sum21')                               8d26s24
          nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv2)                        8d27s24
          do ii=0,nbn(lsb)-1                                            8d26s24
           iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)                   8d26s24
           iadtt=itt4v(lsb,6)+nvirt(isbv1)*iq                           8d26s24
           do jj=0,nan(lsa)-1                                            5d31s24
            jv1=ibc(iptap(lsa)+jj)                                       5d31s24
            do jv2=0,nbasisp(jsbv2)-1                                    6d7s24
             lsum=isum+nvv*(jv2+nbasisp(jsbv2)*jj)                      8d26s24
             iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv2                       8d27s24
     $                    +nbasisp(jsbv2)*noc(isbv2)))                  8d27s24
             do iv2=0,nvirt(isbv2)-1                                     6d6s24
              jadv=iadv+nrowi*iv2                                       8d27s24
              iv1top=iv2+1-ist                                           6d5s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
              irec=nvirt(isbv1)*iv2                                     8d23s24
              itri=((iv2*(iv2+isub))/2)                                 8d23s24
              itri=lsum+itri+isw*(irec-itri)                                 6d5s24
              do iv1=0,iv1top                                            6d5s24
               bc(iadtt+iv1)=bc(iadtt+iv1)+bc(itri+iv1)*eraw2(jadv)     8d23s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          ibcoff=itmp                                                   8d26s24
          call second(time3)
          telapo(6)=telapo(6)+time3-time1-tovr
          end if                                                        8d26s24
         end if                                                         5d31s24
         isub=-1                                                        5d31s24
         nvv=nvvt                                                       5d31s24
        end do                                                          5d31s24
       end if                                                           5d31s24
      end if                                                            5d31s24
      if(.not.lsame)then                                                5d31s24
       isbv1=lsa                                                        5d31s24
       isbv2=lsd                                                        5d31s24
       jsbv1=lsb                                                        5d31s24
       jsbv2=lsc                                                        5d31s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         5d31s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
        if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
         nv1=nvirt(isbv1)                                               6d5s24
         nv2=nvirt(isbv2)                                               6d5s24
         mv1=nbasisp(jsbv1)                                             6d10s24
         mv2=nbasisp(jsbv2)                                             6d10s24
         if(isbv2.eq.isbv1)then                                         6d7s24
          isw=0                                                          5d20s24
          nvvs=(nv1*(nv1+1))/2                                          5d31s24
          nvvt=(nv1*(nv1-1))/2                                          5d31s24
         else                                                            5d20s24
          isw=1                                                          5d20s24
          nvvs=nv1*nv2                                                  5d31s24
          nvvt=nvvs                                                      5d20s24
         end if                                                          5d20s24
         isub=1                                                          5d31s24
         nvv=nvvs                                                        5d31s24
         mvv=mv1*mv2                                                    6d7s24
         nrow=nbn(lsb)*nbasisp(jsbv2)                                   8d26s24
         do ist=1,2                                                      5d31s24
          if(min(nvv,naddr(jtv12,ist)).gt.0)then                        8d27s24
           if(itt4v(lsa,2).gt.0)then                                    8d27s24
           call second(time1)
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv1=ibc(iptbp(lsb)+jj)                                       5d31s24
            do jv2=0,nbasisp(jsbv2)-1                                   6d7s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do iv2=0,nvirt(isbv2)-1                                    6d5s24
              iv2p=iv2+noc(isbv2)                                       6d5s24
              iv1top=iv2+1-ist                                           5d31s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
              do iv1=0,iv1top                                            5d31s24
               irec=iv1+nvirt(isbv1)*iv2                                6d5s24
               itri=((iv2*(iv2+isub))/2)+iv1                             5d31s24
               itri=itri+isw*(irec-itri)                                 5d31s24
               sum=0d0                                                  5d31s24
               iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri             6d11s24
               do i=0,naddr(jtv12,ist)-1                                5d31s24
                sum=sum+bc(jad+i)*bc(iad+i)                             6d11s24
               end do                                                   5d31s24
               do ii=0,nan(lsa)-1                                       5d31s24
                iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)              5d31s24
                iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv2                    5d31s24
     $                    +nbasisp(jsbv2)*iv2p))                           5d21s24
                iadtt=itt4v(lsa,2)+iv1+nvirt(isbv1)*iq                  6d7s24
                bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                     5d31s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           call second(time2)
           telapo(2)=telapo(2)+time2-time1-tovr
           else
           call second(time2)
           end if                                                       8d27s24
           if(itt4v(lsa,6).gt.0)then                                    8d27s24
           itmp=ibcoff                                                  8d26s24
           isum=itmp+nrow*naddr(jtv12,ist)                              8d26s24
           ibcoff=isum+nvv*nrow                                         8d26s24
           call enough('tt4v.sum22',bc,ibc)                             8d23s24
           jtmp=itmp                                                    8d26s24
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv1=ibc(iptbp(lsb)+jj)                                       5d31s24
            do jv2=0,nbasisp(jsbv2)-1                                   6d7s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do i=0,naddr(jtv12,ist)-1                                  8d26s24
              bc(jtmp+i)=bc(jad+i)                                      8d26s24
             end do                                                     8d26s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d26s24
            end do                                                      8d26s24
           end do                                                       8d26s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),1d0,            8d26s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d26s24
     $         bc(isum),nvv,'tt4v.Sum22')                               8d26s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv2)                       8d26s24
           do ii=0,nan(lsa)-1                                           8d26s24
            iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)                  8d26s24
            iadtt=itt4v(lsa,6)+nvirt(isbv1)*iq                          8d26s24
            do jj=0,nbn(lsb)-1                                            5d31s24
             jv1=ibc(iptbp(lsb)+jj)                                       5d31s24
             do jv2=0,nbasisp(jsbv2)-1                                   6d7s24
              lsum=isum+nvv*(jv2+nbasisp(jsbv2)*jj)                     8d26s24
              iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv2                      8d26s24
     $                    +nbasisp(jsbv2)*noc(isbv2)))                  8d26s24
              do iv2=0,nvirt(isbv2)-1                                    6d5s24
               iv1top=iv2+1-ist                                           5d31s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
               irec=nvirt(isbv1)*iv2                                    8d23s24
               itri=((iv2*(iv2+isub))/2)                                8d23s24
               itri=lsum+itri+isw*(irec-itri)                                 5d31s24
               jadv=iadv+nrowi*iv2                                      8d27s24
               do iv1=0,iv1top                                            5d31s24
                bc(iadtt+iv1)=bc(iadtt+iv1)+bc(itri+iv1)*eraw2(jadv)    8d27s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           ibcoff=itmp                                                  8d26s24
           call second(time3)
           telapo(6)=telapo(6)+time3-time2-tovr
           end if                                                       8d26s24
          end if                                                         5d31s24
          isub=-1                                                       5d31s24
          nvv=nvvt                                                      5d31s24
         end do                                                          5d31s24
        end if                                                           5d31s24
       end if                                                            5d31s24
      end if                                                            5d31s24
cc    j3
      isbv2=lsa                                                         5d31s24
      isbv1=lsd                                                         5d31s24
      jsbv2=lsb                                                         5d31s24
      jsbv1=lsc                                                         5d31s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d6s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
       if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
        nv1=nvirt(isbv1)                                                6d5s24
        nv2=nvirt(isbv2)                                                6d5s24
        mv1=nbasisp(jsbv1)                                              6d10s24
        mv2=nbasisp(jsbv2)                                              6d10s24
        if(isbv2.eq.isbv1)then                                          6d7s24
         isw=0                                                          5d20s24
         nvvs=(nv1*(nv1+1))/2                                           5d31s24
         nvvt=(nv1*(nv1-1))/2                                           5d31s24
        else                                                            5d20s24
         isw=1                                                          5d20s24
         nvvs=nv1*nv2                                                   5d31s24
         nvvt=nvvs                                                      5d20s24
        end if                                                          5d20s24
        isub=1                                                          5d31s24
        nvv=nvvs                                                        5d31s24
        mvv=mv1*mv2                                                     6d7s24
        do ist=1,2                                                      5d31s24
         if(min(nvv,naddr(jtv12,ist)).gt.0)then                         8d23s24
          if(itt4v(lsa,3).gt.0)then                                     8d26s24
          call second(time1)
          do jj=0,nbn(lsb)-1                                            5d31s24
           jv2=ibc(iptbp(lsb)+jj)                                       6d7s24
           do jv1=0,nbasisp(jsbv1)-1                                    6d7s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do iv2=0,nvirt(isbv2)-1                                     6d6s24
             iv1top=iv2+1-ist                                           6d5s24
             iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
             do iv1=0,iv1top                                            6d5s24
              iv1p=iv1+noc(isbv1)                                        6d5s24
              irec=iv1+nvirt(isbv1)*iv2                                 6d5s24
              itri=((iv2*(iv2+isub))/2)+iv1                             6d5s24
              itri=itri+isw*(irec-itri)                                 6d5s24
              sum=0d0                                                   5d31s24
              iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri              6d11s24
              do i=0,naddr(jtv12,ist)-1                                 5d31s24
               sum=sum+bc(jad+i)*bc(iad+i)                                  5d31s24
              end do                                                    5d31s24
              do ii=0,nan(lsa)-1                                        5d31s24
               iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)               5d31s24
               iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv1                     6d10s24
     $                    +nbasisp(jsbv1)*iv1p))                        6d10s24
               iadtt=itt4v(lsa,3)+iv2+nvirt(isbv2)*iq                   6d7s24
               bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                      5d31s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          call second(time2)
          telapo(3)=telapo(3)+time2-time1-tovr
          else                                                          8d26s24
           call second(time2)
          end if                                                        8d26s24
          if(itt4v(lsa,7).gt.0)then                                     8d26s24
          itmp=ibcoff                                                   8d23s24
          nrow=nbn(lsb)*nbasisp(jsbv1)                                  8d23s24
          isum=itmp+nrow*naddr(jtv12,ist)                               8d23s24
          ibcoff=isum+nvv*nrow                                          8d23s24
          call enough('tt4v.sum31',bc,ibc)                              8d23s24
          jtmp=itmp                                                     8d23s24
          do jj=0,nbn(lsb)-1                                            5d31s24
           jv2=ibc(iptbp(lsb)+jj)                                       6d7s24
           do jv1=0,nbasisp(jsbv1)-1                                    6d7s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do i=0,naddr(jtv12,ist)-1                                   8d23s24
             bc(jtmp+i)=bc(jad+i)                                       8d27s24
            end do                                                      8d23s24
            jtmp=jtmp+naddr(jtv12,ist)                                  8d27s24
           end do                                                       8d23s24
          end do                                                        8d23s24
          call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),1d0,             8d27s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d27s24
     $         bc(isum),nvv,'tt4v.sum31')                               8d27s24
          nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv1)                        8d23s24
          do jj=0,nbn(lsb)-1                                            5d31s24
           jv2=ibc(iptbp(lsb)+jj)                                       6d7s24
           do ii=0,nan(lsa)-1                                           8d23s24
            iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)                  8d23s24
            iadtt=itt4v(lsa,7)+nvirt(isbv2)*iq                          8d23s24
            do jv1=0,nbasisp(jsbv1)-1                                    6d7s24
             lsum=isum+nvv*(jv1+nbasisp(jsbv1)*jj)                      8d27s24
             iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv1                       8d23s24
     $                    +nbasisp(jsbv1)*noc(isbv1)))                  8d23s24
             do iv2=0,nvirt(isbv2)-1                                     6d6s24
              iv1top=iv2+1-ist                                           6d5s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
              irec=nvirt(isbv1)*iv2                                     8d23s24
              itri=((iv2*(iv2+isub))/2)                                 8d23s24
              itri=lsum+itri+isw*(irec-itri)                            8d27s24
              do iv1=0,iv1top                                            6d5s24
               bc(iadtt+iv2)=bc(iadtt+iv2)                              8d23s24
     $              +bc(itri+iv1)*eraw2(iadv+nrowi*iv1)                 8d27s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          ibcoff=itmp                                                   8d23s24
          call second(time3)
          telapo(7)=telapo(7)+time3-time2-tovr
          end if                                                        8d26s24
         end if                                                         5d31s24
         isub=-1                                                        5d31s24
         nvv=nvvt                                                       5d31s24
        end do                                                          5d31s24
       end if                                                           5d31s24
      end if                                                            5d31s24
      if(.not.lsame)then                                                5d31s24
       isbv2=lsb                                                        5d31s24
       isbv1=lsd                                                        5d31s24
       jsbv2=lsa                                                        5d31s24
       jsbv1=lsc                                                        5d31s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         5d31s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
        if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
         nv1=nvirt(isbv1)                                               6d5s24
         nv2=nvirt(isbv2)                                               6d5s24
         mv1=nbasisp(jsbv1)                                             6d10s24
         mv2=nbasisp(jsbv2)                                             6d10s24
         if(isbv2.eq.isbv1)then                                         6d7s24
          isw=0                                                          5d20s24
          nvvs=(nv1*(nv1+1))/2                                          5d31s24
          nvvt=(nv1*(nv1-1))/2                                          5d31s24
         else                                                            5d20s24
          isw=1                                                          5d20s24
          nvvs=nv1*nv2                                                  5d31s24
          nvvt=nvvs                                                      5d20s24
         end if                                                          5d20s24
         isub=1                                                          5d31s24
         nvv=nvvs                                                        5d31s24
         mvv=mv1*mv2                                                    6d7s24
         do ist=1,2                                                      5d31s24
          if(min(nvv,naddr(jtv12,ist)).gt.0)then                        8d23s24
           if(itt4v(lsb,3).gt.0)then                                    8d26s24
           call second(time1)
           do jj=0,nan(lsa)-1                                            5d31s24
            jv2=ibc(iptap(lsa)+jj)                                      6d7s24
            do jv1=0,nbasisp(jsbv1)-1                                   6d7s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do iv2=0,nvirt(isbv2)-1                                    6d5s24
              iv1top=iv2+1-ist                                           5d31s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
              do iv1=0,iv1top                                            5d31s24
               iv1p=iv1+noc(isbv1)                                       6d5s24
               irec=iv1+nvirt(isbv1)*iv2                                6d5s24
               itri=((iv2*(iv2+isub))/2)+iv1                             5d31s24
               itri=itri+isw*(irec-itri)                                 5d31s24
               sum=0d0                                                  5d31s24
               iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri             6d11s24
               do i=0,naddr(jtv12,ist)-1                                5d31s24
                sum=sum+bc(jad+i)*bc(iad+i)                             6d11s24
               end do                                                   5d31s24
               do ii=0,nbn(lsb)-1                                       5d31s24
                iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)              5d31s24
                iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv1                    6d10s24
     $                    +nbasisp(jsbv1)*iv1p))                        6d10s24
                iadtt=itt4v(lsb,3)+iv2+nvirt(isbv2)*iq                  6d7s24
                bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                     5d31s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           call second(time2)
           telapo(3)=telapo(3)+time2-time1-tovr
           else                                                         8d26s24
            call second(time2)
           end if                                                       8d26s24
           if(itt4v(lsb,7).gt.0)then                                    8d26s24
           itmp=ibcoff                                                  8d23s24
           nrow=nan(lsa)*nbasisp(jsbv1)                                 8d23s24
           isum=itmp+nrow*naddr(jtv12,ist)                              8d23s24
           ibcoff=isum+nrow*nvv                                         8d23s24
           call enough('tt4v.sum32',bc,ibc)                             8d23s24
           jtmp=itmp                                                    8d23s24
           do jj=0,nan(lsa)-1                                            5d31s24
            jv2=ibc(iptap(lsa)+jj)                                      6d7s24
            do jv1=0,nbasisp(jsbv1)-1                                   6d7s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do i=0,naddr(jtv12,ist)-1                                  8d23s24
              bc(jtmp+i)=bc(jad+i)                                      8d27s24
             end do                                                     8d23s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d27s24
            end do                                                      8d23s24
           end do                                                       8d23s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),1d0,            8d27s24
     $        bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0, 8d27s24
     $        bc(isum),nvv,'tt4v.sum32')                                8d27s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv1)                       8d23s24
           do jj=0,nan(lsa)-1                                            5d31s24
            jv2=ibc(iptap(lsa)+jj)                                      6d7s24
            do ii=0,nbn(lsb)-1                                          8d23s24
             iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)                 8d23s24
             iadtt=itt4v(lsb,7)+nvirt(isbv2)*iq                         8d23s24
             do jv1=0,nbasisp(jsbv1)-1                                   6d7s24
              lsum=isum+nvv*(jv1+nbasisp(jsbv1)*jj)                     8d27s24
              iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv1                      8d23s24
     $                    +nbasisp(jsbv1)*noc(isbv1)))                  8d23s24
              do iv2=0,nvirt(isbv2)-1                                    6d5s24
               iv1top=iv2+1-ist                                         8d23s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                8d23s24
               irec=nvirt(isbv1)*iv2                                    8d23s24
               itri=((iv2*(iv2+isub))/2)                                8d23s24
               itri=lsum+itri+isw*(irec-itri)                           8d27s24
               do iv1=0,iv1top                                            5d31s24
                bc(iadtt+iv2)=bc(iadtt+iv2)                             8d23s24
     $               +bc(itri+iv1)*eraw2(iadv+nrowi*iv1)                8d27s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           ibcoff=itmp                                                  8d23s24
           call second(time3)
           telapo(7)=telapo(7)+time3-time2-tovr
           end if                                                       8d26s24
          end if                                                         5d31s24
          isub=-1                                                       5d31s24
          nvv=nvvt                                                      5d31s24
         end do                                                          5d31s24
        end if                                                           5d31s24
       end if                                                            5d31s24
      end if                                                            5d31s24
cc    j4
      isbv2=lsb                                                         6d10s24
      isbv1=lsd                                                         6d10s24
      jsbv2=lsa                                                         6d10s24
      jsbv1=lsc                                                         6d10s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d6s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
       if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
        nv1=nvirt(isbv1)                                                6d5s24
        nv2=nvirt(isbv2)                                                6d5s24
        mv1=nbasisp(jsbv1)                                              6d10s24
        mv2=nbasisp(jsbv2)                                              6d10s24
        if(isbv2.eq.isbv1)then                                          6d7s24
         isw=0                                                          5d20s24
         nvvs=(nv1*(nv1+1))/2                                           5d31s24
         nvvt=(nv1*(nv1-1))/2                                           5d31s24
        else                                                            5d20s24
         isw=1                                                          5d20s24
         nvvs=nv1*nv2                                                   5d31s24
         nvvt=nvvs                                                      5d20s24
        end if                                                          5d20s24
        isub=1                                                          5d31s24
        nvv=nvvs                                                        5d31s24
        mvv=mv1*mv2                                                     6d7s24
        do ist=1,2                                                      5d31s24
         if(min(nvv,naddr(jtv12,ist)).gt.0)then                         8d23s24
          if(itt4v(lsb,4).gt.0)then                                     8d26s24
          call second(time1)
          do jj=0,nan(lsa)-1                                            5d31s24
           jv2=ibc(iptap(lsa)+jj)                                       6d10s24
           do jv1=0,nbasisp(jsbv1)-1                                    6d10s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do iv2=0,nvirt(isbv2)-1                                     6d6s24
             iv1top=iv2+1-ist                                           6d5s24
             iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
             do iv1=0,iv1top                                            6d5s24
              iv1p=iv1+noc(isbv1)                                       6d10s24
              irec=iv1+nvirt(isbv1)*iv2                                 6d5s24
              itri=((iv2*(iv2+isub))/2)+iv1                             6d5s24
              itri=itri+isw*(irec-itri)                                 6d5s24
              sum=0d0                                                   5d31s24
              iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri              6d11s24
              do i=0,naddr(jtv12,ist)-1                                 5d31s24
               sum=sum+bc(jad+i)*bc(iad+i)                              6d11s24
              end do                                                    5d31s24
              do ii=0,nbn(lsb)-1                                        5d31s24
               iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)               5d31s24
               iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv1                     6d10s24
     $                    +nbasisp(jsbv1)*iv1p))                        6d10s24
               iadtt=itt4v(lsb,4)+iv2+nvirt(isbv2)*iq                   6d10s24
               bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                      5d31s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          call second(time2)
          telapo(4)=telapo(4)+time2-time1-tovr
          else                                                          8d26s24
          call second(time2)
          end if                                                        8d26s24
          if(itt4v(lsb,8).gt.0)then                                     8d26s24
          nrow=nan(lsa)*nbasisp(jsbv1)                                  8d23s24
          itmp=ibcoff                                                   8d23s24
          isum=itmp+nrow*naddr(jtv12,ist)                               8d23s24
          ibcoff=isum+nvv*nrow                                          8d23s24
          call enough('tt4v.sum41',bc,ibc)                              8d23s24
          jtmp=itmp                                                     8d23s24
          do jj=0,nan(lsa)-1                                            5d31s24
           jv2=ibc(iptap(lsa)+jj)                                       6d10s24
           do jv1=0,nbasisp(jsbv1)-1                                    6d10s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do i=0,naddr(jtv12,ist)-1                                   8d23s24
             bc(jtmp+i)=bc(jad+i)                                       8d23s24
            end do                                                      8d23s24
            jtmp=jtmp+naddr(jtv12,ist)                                  8d23s24
           end do                                                       8d23s24
          end do                                                        8d23s24
          call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),1d0,             8d23s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d23s24
     $         bc(isum),nvv,'tt4v.sum41')                               8d23s24
          nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv1)                        8d23s24
          do jj=0,nan(lsa)-1                                            5d31s24
           jv2=ibc(iptap(lsa)+jj)                                       6d10s24
           do ii=0,nbn(lsb)-1                                           8d23s24
            iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)                  8d23s24
            iadtt=itt4v(lsb,8)+nvirt(isbv2)*iq                          8d23s24
            do jv1=0,nbasisp(jsbv1)-1                                    6d10s24
             lsum=isum+nvv*(jv1+nbasisp(jsbv1)*jj)                      8d23s24
             iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv1                       8d23s24
     $                    +nbasisp(jsbv1)*noc(isbv1)))                        6d10s24
             do iv2=0,nvirt(isbv2)-1                                     6d6s24
              iv1top=iv2+1-ist                                           6d5s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
              irec=nvirt(isbv1)*iv2                                     8d23s24
              itri=((iv2*(iv2+isub))/2)                                 8d23s24
              itri=lsum+itri+isw*(irec-itri)                            8d23s24
              do iv1=0,iv1top                                            6d5s24
               bc(iadtt+iv2)=bc(iadtt+iv2)                              8d23s24
     $              +bc(itri+iv1)*eraw2(iadv+nrowi*iv1)                 8d23s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          ibcoff=itmp                                                   8d23s24
          call second(time3)
          telapo(8)=telapo(8)+time3-time2-tovr
          end if                                                        8d26s24
         end if                                                         5d31s24
         isub=-1                                                        5d31s24
         nvv=nvvt                                                       5d31s24
        end do                                                          5d31s24
       end if                                                           5d31s24
      end if                                                            5d31s24
      if(.not.lsame)then                                                5d31s24
       isbv2=lsa                                                        6d10s24
       isbv1=lsd                                                        6d10s24
       jsbv2=lsb                                                        6d10s24
       jsbv1=lsc                                                        6d10s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         5d31s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
        if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
         nv1=nvirt(isbv1)                                               6d5s24
         nv2=nvirt(isbv2)                                               6d5s24
         mv1=nbasisp(jsbv1)                                             6d10s24
         mv2=nbasisp(jsbv2)                                             6d10s24
         if(isbv2.eq.isbv1)then                                         6d7s24
          isw=0                                                          5d20s24
          nvvs=(nv1*(nv1+1))/2                                          5d31s24
          nvvt=(nv1*(nv1-1))/2                                          5d31s24
         else                                                            5d20s24
          isw=1                                                          5d20s24
          nvvs=nv1*nv2                                                  5d31s24
          nvvt=nvvs                                                      5d20s24
         end if                                                          5d20s24
         isub=1                                                          5d31s24
         nvv=nvvs                                                        5d31s24
         mvv=mv1*mv2                                                    6d7s24
         do ist=1,2                                                      5d31s24
          if(min(nvv,naddr(jtv12,ist)).gt.0)then                        8d23s24
           if(itt4v(lsa,4).gt.0)then                                    8d26s24
           call second(time1)
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv2=ibc(iptbp(lsb)+jj)                                      6d10s24
            do jv1=0,nbasisp(jsbv1)-1                                   6d10s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do iv2=0,nvirt(isbv2)-1                                    6d5s24
              iv1top=iv2+1-ist                                           5d31s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
              do iv1=0,iv1top                                            5d31s24
               iv1p=iv1+noc(isbv1)                                      6d10s24
               irec=iv1+nvirt(isbv1)*iv2                                6d5s24
               itri=((iv2*(iv2+isub))/2)+iv1                             5d31s24
               itri=itri+isw*(irec-itri)                                 5d31s24
               sum=0d0                                                  5d31s24
               iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri             6d11s24
               do i=0,naddr(jtv12,ist)-1                                5d31s24
                sum=sum+bc(jad+i)*bc(iad+i)                             6d11s24
               end do                                                   5d31s24
               do ii=0,nan(lsa)-1                                       5d31s24
                iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)              5d31s24
                iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv1                    6d10s24
     $                    +nbasisp(jsbv1)*iv1p))                        6d10s24
                iadtt=itt4v(lsa,4)+iv2+nvirt(isbv2)*iq                  6d10s24
                bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                     5d31s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           call second(time2)
           telapo(4)=telapo(4)+time2-time1-tovr
           else                                                         8d26s24
            call second(time2)
           end if                                                       8d26s24
           if(itt4v(lsa,8).gt.0)then                                    8d26s24
           nrow=nbn(lsb)*nbasisp(jsbv1)                                 8d23s24
           itmp=ibcoff                                                  8d23s24
           isum=itmp+nrow*naddr(jtv12,ist)                              8d23s24
           ibcoff=isum+nvv*nrow                                         8d23s24
           call enough('tt4v.sum42',bc,ibc)                             8d23s24
           jtmp=itmp                                                    8d23s24
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv2=ibc(iptbp(lsb)+jj)                                      6d10s24
            do jv1=0,nbasisp(jsbv1)-1                                   6d10s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do i=0,naddr(jtv12,ist)-1                                  8d23s24
              bc(jtmp+i)=bc(jad+i)                                      8d23s24
             end do                                                     8d23s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d23s24
            end do                                                      8d23s24
           end do                                                       8d23s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),1d0,            8d23s24
     $       bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,  8d23s24
     $          bc(isum),nvv,'tt4v.sum42')                              8d23s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv1)                       8d23s24
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv2=ibc(iptbp(lsb)+jj)                                      6d10s24
            do ii=0,nan(lsa)-1                                          8d23s24
             iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)                 8d23s24
             iadtt=itt4v(lsa,8)+nvirt(isbv2)*iq                         8d23s24
             do jv1=0,nbasisp(jsbv1)-1                                   6d10s24
              lsum=isum+nvv*(jv1+nbasisp(jsbv1)*jj)                     8d23s24
              iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv1                      8d23s24
     $                    +nbasisp(jsbv1)*noc(isbv1)))                        6d10s24
              do iv2=0,nvirt(isbv2)-1                                    6d5s24
               iv1top=iv2+1-ist                                           5d31s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
               irec=nvirt(isbv1)*iv2                                    8d23s24
               itri=((iv2*(iv2+isub))/2)                                8d23s24
               itri=lsum+itri+isw*(irec-itri)                           8d23s24
               do iv1=0,iv1top                                            5d31s24
                bc(iadtt+iv2)=bc(iadtt+iv2)                             8d23s24
     $               +bc(itri+iv1)*eraw2(iadv+nrowi*iv1)                8d23s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           ibcoff=itmp                                                  8d23s24
           call second(time3)
           telapo(8)=telapo(8)+time3-time2-tovr
           end if                                                       8d26s24
          end if                                                         5d31s24
          isub=-1                                                       5d31s24
          nvv=nvvt                                                      5d31s24
         end do                                                          5d31s24
        end if                                                           5d31s24
       end if                                                            5d31s24
      end if                                                            5d31s24
cc k 1
      isbv1=lsa                                                         5d31s24
      isbv2=lsd                                                         5d31s24
      jsbv2=lsb                                                         5d31s24
      jsbv1=lsc                                                         5d31s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d6s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
       if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
        nv1=nvirt(isbv1)                                                6d5s24
        nv2=nvirt(isbv2)                                                6d5s24
        mv1=nbasisp(jsbv1)                                              6d10s24
        mv2=nbasisp(jsbv2)                                              6d10s24
        if(isbv2.eq.isbv1)then                                          6d7s24
         isw=0                                                          5d20s24
         nvvs=(nv1*(nv1+1))/2                                           5d31s24
         nvvt=(nv1*(nv1-1))/2                                           5d31s24
        else                                                            5d20s24
         isw=1                                                          5d20s24
         nvvs=nv1*nv2                                                   5d31s24
         nvvt=nvvs                                                      5d20s24
        end if                                                          5d20s24
        isub=1                                                          5d31s24
        nvv=nvvs                                                        5d31s24
        mvv=mv1*mv2                                                     6d7s24
        phase=1d0                                                       6d11s24
        nrow=nbn(lsb)*nbasisp(jsbv1)                                    8d26s24
        do ist=1,2                                                      5d31s24
         if(min(nvv,naddr(jtv12,ist)).gt.0)then                         8d26s24
          if(itt4v(lsa,1).gt.0)then                                     8d26s24
          call second(time1)
          do jj=0,nbn(lsb)-1                                            5d31s24
           jv2=ibc(iptbp(lsb)+jj)                                       6d11s24
           do jv1=0,nbasisp(jsbv1)-1                                    6d11s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do iv2=0,nvirt(isbv2)-1                                     6d6s24
             iv1top=iv2+1-ist                                           6d5s24
             iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
             iv2p=iv2+noc(isbv2)                                        6d5s24
             do iv1=0,iv1top                                            6d5s24
              irec=iv1+nvirt(isbv1)*iv2                                 6d5s24
              itri=((iv2*(iv2+isub))/2)+iv1                             6d5s24
              itri=itri+isw*(irec-itri)                                 6d5s24
              iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri              6d11s24
              sum=0d0                                                   6d11s24
              do i=0,naddr(jtv12,ist)-1                                 6d11s24
               orig=sum
               sum=sum+bc(jad+i)*bc(iad+i)                              6d11s24
              end do                                                    6d11s24
              sum=sum*phase                                             6d11s24
              do ii=0,nan(lsa)-1                                        5d31s24
               iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)               5d31s24
               iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv1                     6d11s24
     $                    +nbasisp(jsbv1)*iv2p))                        6d11s24
               iadtt=itt4v(lsa,1)+iv1+nvirt(isbv1)*iq                   6d7s24
               bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                      5d31s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          call second(time2)
          telapo(1)=telapo(1)+time2-time1-tovr
          else                                                          8d26s24
           call second(time2)
          end if                                                        8d26s24
          if(itt4v(lsa,5).gt.0)then                                     8d27s24
           itmp=ibcoff                                                  8d26s24
           isum=itmp+nrow*naddr(jtv12,ist)                              8d26s24
           ibcoff=isum+nrow*nvv                                         8d26s24
           call enough('tt4v.sum13b',bc,ibc)                            8d26s24
           jtmp=itmp                                                    8d26s24
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv2=ibc(iptbp(lsb)+jj)                                       6d11s24
            do jv1=0,nbasisp(jsbv1)-1                                    6d11s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
             do i=0,naddr(jtv12,ist)-1                                  8d26s24
              bc(jtmp+i)=bc(jad+i)                                      8d26s24
             end do                                                     8d26s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d26s24
            end do                                                      8d26s24
           end do                                                       8d26s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),phase,          8d26s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d26s24
     $         bc(isum),nvv,'tt4v.sum13b')                              8d26s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv1)                       8d26s24
           do ii=0,nan(lsa)-1                                           8d26s24
            iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)                  8d26s24
            iadtt=itt4v(lsa,5)+nvirt(isbv1)*iq                          8d26s24
            do jj=0,nbn(lsb)-1                                            5d31s24
             jv2=ibc(iptbp(lsb)+jj)                                       6d11s24
             do jv1=0,nbasisp(jsbv1)-1                                    6d11s24
              lsum=isum+nvv*(jv1+nbasisp(jsbv1)*jj)                     8d26s24
              iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv1                      8d23s24
     $                    +nbasisp(jsbv1)*noc(isbv2)))                        6d11s24
              do iv2=0,nvirt(isbv2)-1                                     6d6s24
               iv1top=iv2+1-ist                                           6d5s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                8d27s24
               irec=nvirt(isbv1)*iv2                                    8d27s24
               itri=((iv2*(iv2+isub))/2)                                8d27s24
               itri=lsum+itri+isw*(irec-itri)                           8d27s24
               jadv=iadv+nrowi*iv2
               do iv1=0,iv1top                                           8d23s24
                bc(iadtt+iv1)=bc(iadtt+iv1)+bc(itri+iv1)*eraw2(jadv)    8d26s24
               end do                                                   8d26s24
              end do                                                    8d26s24
             end do                                                     8d26s24
            end do                                                      8d26s24
           end do                                                       8d26s24
           ibcoff=itmp                                                  8d26s24
           call second(time3)
           telapo(5)=telapo(5)+time3-time2-tovr
          end if                                                        8d26s24
         end if                                                         5d31s24
         isub=-1                                                        5d31s24
         nvv=nvvt                                                       5d31s24
         phase=-1d0                                                     6d11s24
        end do                                                          5d31s24
       end if                                                           5d31s24
      end if                                                            5d31s24
      if(.not.lsame)then                                                5d31s24
       isbv1=lsb                                                        5d31s24
       isbv2=lsd                                                        5d31s24
       jsbv2=lsa                                                        5d31s24
       jsbv1=lsc                                                        5d31s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         5d31s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
        if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
         nv1=nvirt(isbv1)                                               6d5s24
         nv2=nvirt(isbv2)                                               6d5s24
         mv1=nbasisp(jsbv1)                                             6d10s24
         mv2=nbasisp(jsbv2)                                             6d10s24
         if(isbv2.eq.isbv1)then                                         6d7s24
          isw=0                                                          5d20s24
          nvvs=(nv1*(nv1+1))/2                                          5d31s24
          nvvt=(nv1*(nv1-1))/2                                          5d31s24
         else                                                            5d20s24
          isw=1                                                          5d20s24
          nvvs=nv1*nv2                                                  5d31s24
          nvvt=nvvs                                                      5d20s24
         end if                                                          5d20s24
         isub=1                                                          5d31s24
         nvv=nvvs                                                        5d31s24
         mvv=mv1*mv2                                                    6d7s24
         phase=1d0                                                      6d11s24
         nrow=nan(lsa)*nbasisp(jsbv1)                                   8d26s24
         do ist=1,2                                                      5d31s24
          if(min(nvv,naddr(jtv12,ist)).gt.0)then                        8d26s24
           if(itt4v(lsb,1).gt.0)then                                    8d26s24
           call second(time1)
           do jj=0,nan(lsa)-1                                            5d31s24
            jv2=ibc(iptap(lsa)+jj)                                      6d11s24
            do jv1=0,nbasisp(jsbv1)-1                                   6d11s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do iv2=0,nvirt(isbv2)-1                                    6d5s24
              iv2p=iv2+noc(isbv2)                                       6d5s24
              iv1top=iv2+1-ist                                           5d31s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
              do iv1=0,iv1top                                            5d31s24
               irec=iv1+nvirt(isbv1)*iv2                                6d5s24
               itri=((iv2*(iv2+isub))/2)+iv1                             5d31s24
               itri=itri+isw*(irec-itri)                                 5d31s24
               sum=0d0                                                  5d31s24
               iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri             6d11s24
               do i=0,naddr(jtv12,ist)-1                                5d31s24
                sum=sum+bc(jad+i)*bc(iad+i)                             6d11s24
               end do                                                   5d31s24
               sum=sum*phase                                            6d11s24
               do ii=0,nbn(lsb)-1                                       5d31s24
                iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)              5d31s24
                iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv1                    6d11s24
     $                    +nbasisp(jsbv1)*iv2p))                        6d11s24
                iadtt=itt4v(lsb,1)+iv1+nvirt(isbv1)*iq                  6d7s24
                bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                     5d31s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           call second(time2)
           telapo(1)=telapo(1)+time2-time1-tovr
           else                                                         8d26s24
            call second(time2)
           end if                                                       8d26s24
           if(itt4v(lsb,5).gt.0)then                                    8d27s24
            itmp=ibcoff                                                 8d26s24
            isum=itmp+nrow*naddr(jtv12,ist)                             8d26s24
            ibcoff=isum+nvv*nrow                                        8d26s24
            call enough('tt4v.sum14',bc,ibc)                            8d26s24
            jtmp=itmp                                                   8d26s24
           do jj=0,nan(lsa)-1                                            5d31s24
            jv2=ibc(iptap(lsa)+jj)                                      6d11s24
            do jv1=0,nbasisp(jsbv1)-1                                   6d11s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            8d26s24
             do i=0,naddr(jtv12,ist)-1                                  8d26s24
              bc(jtmp+i)=bc(jad+i)                                      8d26s24
             end do                                                     8d26s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d26s24
            end do                                                      8d26s24
           end do                                                       8d26s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),phase,          8d26s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d26s24
     $          bc(isum),nvv,'tt4v.sum14b')                              8d26s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv1)                       8d26s24
           do ii=0,nbn(lsb)-1                                           8d26s24
            iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)                  8d26s24
            iadtt=itt4v(lsb,5)+nvirt(isbv1)*iq                          8d26s24
            do jj=0,nan(lsa)-1                                            5d31s24
             jv2=ibc(iptap(lsa)+jj)                                      6d11s24
             do jv1=0,nbasisp(jsbv1)-1                                   6d11s24
              iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv1                      8d26s24
     $                    +nbasisp(jsbv1)*noc(isbv2)))                        6d11s24
              lsum=isum+nvv*(jv1+nbasisp(jsbv1)*jj)                     8d26s24
              do iv2=0,nvirt(isbv2)-1                                    6d5s24
               jadv=iadv+nrowi*iv2                                      8d26s24
               irec=nvirt(isbv1)*iv2                                    8d27s24
               itri=((iv2*(iv2+isub))/2)                                8d27s24
               itri=lsum+itri+isw*(irec-itri)                           8d27s24
               iv1top=iv2+1-ist                                           5d31s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                8d27s24
               do iv1=0,iv1top                                          8d26s24
                bc(iadtt+iv1)=bc(iadtt+iv1)+bc(itri+iv1)*eraw2(jadv)    8d23s24
               end do                                                   8d26s24
              end do                                                    8d26s24
             end do                                                     8d26s24
            end do                                                      8d26s24
           end do                                                       8d26s24
           ibcoff=itmp                                                  8d26s24
           call second(time3)
           telapo(5)=telapo(5)+time3-time2-tovr
           end if                                                       8d26s24
          end if                                                         5d31s24
          isub=-1                                                       5d31s24
          nvv=nvvt                                                      5d31s24
          phase=-1d0                                                    6d11s24
         end do                                                          5d31s24
        end if                                                           5d31s24
       end if                                                            5d31s24
      end if                                                            5d31s24
cc k2
      isbv1=lsb                                                         5d31s24
      isbv2=lsd                                                         5d31s24
      jsbv2=lsa                                                         5d31s24
      jsbv1=lsc                                                         5d31s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d6s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
       if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
        nv1=nvirt(isbv1)                                                6d5s24
        nv2=nvirt(isbv2)                                                6d5s24
        mv1=nbasisp(jsbv1)                                              6d10s24
        mv2=nbasisp(jsbv2)                                              6d10s24
        if(isbv2.eq.isbv1)then                                          6d7s24
         isw=0                                                          5d20s24
         nvvs=(nv1*(nv1+1))/2                                           5d31s24
         nvvt=(nv1*(nv1-1))/2                                           5d31s24
        else                                                            5d20s24
         isw=1                                                          5d20s24
         nvvs=nv1*nv2                                                   5d31s24
         nvvt=nvvs                                                      5d20s24
        end if                                                          5d20s24
        isub=1                                                          5d31s24
        nvv=nvvs                                                        5d31s24
        mvv=mv1*mv2                                                     6d7s24
        phase=1d0                                                       6d11s24
        nrow=nan(lsa)*nbasisp(jsbv1)                                    8d26s24
        do ist=1,2                                                      5d31s24
         if(min(nvv,naddr(jtv12,ist)).gt.0)then                         8d27s24
          if(itt4v(lsb,2).gt.0)then                                     8d26s24
          call second(time1)
          do jj=0,nan(lsa)-1                                            5d31s24
           jv2=ibc(iptap(lsa)+jj)                                       6d11s24
           do jv1=0,nbasisp(jsbv1)-1                                    6d11s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do iv2=0,nvirt(isbv2)-1                                     6d6s24
             iv1top=iv2+1-ist                                           6d5s24
             iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
             iv2p=iv2+noc(isbv2)                                        6d5s24
             do iv1=0,iv1top                                            6d5s24
              irec=iv1+nvirt(isbv1)*iv2                                 6d5s24
              itri=((iv2*(iv2+isub))/2)+iv1                             6d5s24
              itri=itri+isw*(irec-itri)                                 6d5s24
              sum=0d0                                                   5d31s24
              iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri              6d11s24
              do i=0,naddr(jtv12,ist)-1                                 5d31s24
               sum=sum+bc(jad+i)*bc(iad+i)                              6d11s24
              end do                                                    5d31s24
              sum=sum*phase                                             6d11s24
              do ii=0,nbn(lsb)-1                                        5d31s24
               iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)               5d31s24
               iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv1                     6d11s24
     $                    +nbasisp(jsbv1)*iv2p))                        6d11s24
               iadtt=itt4v(lsb,2)+iv1+nvirt(isbv1)*iq                   6d7s24
               bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                      5d31s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          call second(time2)
          telapo(2)=telapo(2)+time2-time1-tovr
          else
          call second(time2)
          end if                                                        8d26s24
          if(itt4v(lsb,6).gt.0)then                                     8d27s24
           itmp=ibcoff                                                  8d27s24
           isum=itmp+nrow*naddr(jtv12,ist)                              8d27s24
           ibcoff=isum+nvv*nrow                                         8d27s24
           call enough('tt4v.sum23',bc,ibc)                              8d23s24
           jtmp=itmp                                                    8d27s24
           do jj=0,nan(lsa)-1                                            5d31s24
            jv2=ibc(iptap(lsa)+jj)                                       6d11s24
            do jv1=0,nbasisp(jsbv1)-1                                    6d11s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
             do i=0,naddr(jtv12,ist)-1                                  8d27s24
              bc(jtmp+i)=bc(jad+i)                                      8d27s24
             end do                                                     8d27s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d27s24
            end do                                                      8d27s24
           end do                                                       8d27s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),phase,          8d27s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d27s24
     $         bc(isum),nvv,'tt4v.sum23')                               8d27s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv1)                       8d27s24
           do ii=0,nbn(lsb)-1                                           8d23s24
            iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)                  8d23s24
            iadtt=itt4v(lsb,6)+nvirt(isbv1)*iq                          8d23s24
            do jj=0,nan(lsa)-1                                            5d31s24
             do jv1=0,nbasisp(jsbv1)-1                                    6d11s24
              lsum=isum+nvv*(jv1+nbasisp(jsbv1)*jj)                     8d27s24
              iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv1                      8d23s24
     $                    +nbasisp(jsbv1)*noc(isbv2)))                        6d11s24
              do iv2=0,nvirt(isbv2)-1                                     6d6s24
               iv1top=iv2+1-ist                                           6d5s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
               jadv=iadv+nrowi*iv2                                      8d27s24
               irec=nvirt(isbv1)*iv2                                     8d23s24
               itri=((iv2*(iv2+isub))/2)                                 8d23s24
               itri=lsum+itri+isw*(irec-itri)                           8d27s24
               do iv1=0,iv1top                                            6d5s24
                bc(iadtt+iv1)=bc(iadtt+iv1)+bc(itri+iv1)*eraw2(jadv)     8d23s24
               end do                                                    5d31s24
              end do                                                    8d27s24
             end do                                                     8d27s24
            end do                                                      8d27s24
           end do                                                       8d27s24
           ibcoff=itmp                                                  8d27s24
           call second(time3)
           telapo(6)=telapo(6)+time3-time2-tovr
          end if                                                        8d26s24
         end if                                                         5d31s24
         isub=-1                                                        5d31s24
         nvv=nvvt                                                       5d31s24
         phase=-1d0                                                     6d11s24
        end do                                                          5d31s24
       end if                                                           5d31s24
      end if                                                            5d31s24
      if(.not.lsame)then                                                5d31s24
       isbv1=lsa                                                        5d31s24
       isbv2=lsd                                                        5d31s24
       jsbv2=lsb                                                        6d11s24
       jsbv1=lsc                                                        6d11s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         5d31s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
        if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
         nv1=nvirt(isbv1)                                               6d5s24
         nv2=nvirt(isbv2)                                               6d5s24
         mv1=nbasisp(jsbv1)                                             6d10s24
         mv2=nbasisp(jsbv2)                                             6d10s24
         if(isbv2.eq.isbv1)then                                         6d7s24
          isw=0                                                          5d20s24
          nvvs=(nv1*(nv1+1))/2                                          5d31s24
          nvvt=(nv1*(nv1-1))/2                                          5d31s24
         else                                                            5d20s24
          isw=1                                                          5d20s24
          nvvs=nv1*nv2                                                  5d31s24
          nvvt=nvvs                                                      5d20s24
         end if                                                          5d20s24
         isub=1                                                          5d31s24
         nvv=nvvs                                                        5d31s24
         mvv=mv1*mv2                                                    6d7s24
         phase=1d0                                                      6d11s24
         nrow=nbn(lsb)*nbasisp(jsbv1)                                   8d26s24
         do ist=1,2                                                      5d31s24
          if(min(nvv,naddr(jtv12,ist)).gt.0)then                        8d27s24
           if(itt4v(lsa,2).gt.0)then                                    8d26s24
           call second(time1)
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv2=ibc(iptbp(lsb)+jj)                                      6d11s24
            do jv1=0,nbasisp(jsbv1)-1                                   6d11s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do iv2=0,nvirt(isbv2)-1                                    6d5s24
              iv2p=iv2+noc(isbv2)                                       6d5s24
              iv1top=iv2+1-ist                                           5d31s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
              do iv1=0,iv1top                                            5d31s24
               irec=iv1+nvirt(isbv1)*iv2                                6d5s24
               itri=((iv2*(iv2+isub))/2)+iv1                             5d31s24
               itri=itri+isw*(irec-itri)                                 5d31s24
               sum=0d0                                                  5d31s24
               iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri             6d11s24
               do i=0,naddr(jtv12,ist)-1                                5d31s24
                sum=sum+bc(jad+i)*bc(iad+i)                             6d11s24
               end do                                                   5d31s24
               sum=sum*phase                                            6d11s24
               do ii=0,nan(lsa)-1                                       5d31s24
                iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)              5d31s24
                iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv1                    6d11s24
     $                    +nbasisp(jsbv1)*iv2p))                        6d11s24
                iadtt=itt4v(lsa,2)+iv1+nvirt(isbv1)*iq                  6d7s24
                bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                     5d31s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           call second(time2)
           telapo(2)=telapo(2)+time2-time1-tovr
           else                                                         8d27s24
           call second(time2)
           end if                                                       8d26s24
           if(itt4v(lsa,6).gt.0)then                                    8d27s24
            itmp=ibcoff                                                 8d27s24
            isum=itmp+nrow*naddr(jtv12,ist)                             8d27s24
            ibcoff=isum+nrow*nvv                                        8d27s24
            call enough('tt4v.sum24',bc,ibc)                            8d27s24
            jtmp=itmp                                                   8d27s24
            do jj=0,nbn(lsb)-1                                            5d31s24
             jv2=ibc(iptbp(lsb)+jj)                                      6d11s24
             do jv1=0,nbasisp(jsbv1)-1                                   6d11s24
              jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
              jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
              do i=0,naddr(jtv12,ist)-1                                 8d27s24
               bc(jtmp+i)=bc(jad+i)                                     8d27s24
              end do                                                    8d27s24
              jtmp=jtmp+naddr(jtv12,ist)                                8d27s24
             end do                                                     8d27s24
            end do                                                      8d27s24
            call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),phase,         8d27s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d27s24
     $         bc(isum),nvv,'tt4v.sum24')                               8d27s24
            nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv1)                      8d27s24
            do ii=0,nan(lsa)-1                                          8d23s24
             iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)                 8d23s24
             iadtt=itt4v(lsa,6)+nvirt(isbv1)*iq                         8d23s24
             do jj=0,nbn(lsb)-1                                            5d31s24
              do jv1=0,nbasisp(jsbv1)-1                                   6d11s24
               lsum=isum+nvv*(jv1+nbasisp(jsbv1)*jj)                    8d27s24
               iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv1                     8d23s24
     $                    +nbasisp(jsbv1)*noc(isbv2)))                  8d27s24
               do iv2=0,nvirt(isbv2)-1                                   8d23s24
                jadv=iadv+nrowi*iv2                                     8d27s24
                iv1top=iv2+1-ist                                         8d23s24
                iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                8d23s24
                irec=nvirt(isbv1)*iv2                                    8d23s24
                itri=((iv2*(iv2+isub))/2)                                8d23s24
                itri=lsum+itri+isw*(irec-itri)                           8d23s24
                do iv1=0,iv1top                                          8d23s24
                 bc(iadtt+iv1)=bc(iadtt+iv1)+bc(itri+iv1)*eraw2(jadv)    8d23s24
                end do                                                   5d31s24
               end do                                                   8d27s24
              end do                                                    8d27s24
             end do                                                     8d27s24
            end do                                                      8d27s24
            ibcoff=itmp                                                 8d27s24
            call second(time3)
            telapo(6)=telapo(6)+time3-time2-tovr
           end if
          end if                                                         5d31s24
          isub=-1                                                       5d31s24
          nvv=nvvt                                                      5d31s24
          phase=-1d0                                                    6d11s24
         end do                                                          5d31s24
        end if                                                           5d31s24
       end if                                                            5d31s24
      end if                                                            5d31s24
cc k3
      isbv2=lsa                                                         5d31s24
      isbv1=lsd                                                         5d31s24
      jsbv1=lsb                                                         6d11s24
      jsbv2=lsc                                                         6d11s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d6s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
       if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
        nv1=nvirt(isbv1)                                                6d5s24
        nv2=nvirt(isbv2)                                                6d5s24
        mv1=nbasisp(jsbv1)                                              6d10s24
        mv2=nbasisp(jsbv2)                                              6d10s24
        if(isbv2.eq.isbv1)then                                          6d7s24
         isw=0                                                          5d20s24
         nvvs=(nv1*(nv1+1))/2                                           5d31s24
         nvvt=(nv1*(nv1-1))/2                                           5d31s24
        else                                                            5d20s24
         isw=1                                                          5d20s24
         nvvs=nv1*nv2                                                   5d31s24
         nvvt=nvvs                                                      5d20s24
        end if                                                          5d20s24
        isub=1                                                          5d31s24
        nvv=nvvs                                                        5d31s24
        mvv=mv1*mv2                                                     6d7s24
        phase=1d0                                                       6d11s24
        do ist=1,2                                                      5d31s24
         if(min(nvv,naddr(jtv12,ist)).gt.0)then                         8d23s24
          if(itt4v(lsa,3).gt.0)then                                     8d26s24
          call second(time1)
          do jj=0,nbn(lsb)-1                                            5d31s24
           jv1=ibc(iptbp(lsb)+jj)                                       6d11s24
           do jv2=0,nbasisp(jsbv2)-1                                    6d11s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do iv2=0,nvirt(isbv2)-1                                     6d6s24
             iv1top=iv2+1-ist                                           6d5s24
             iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
             do iv1=0,iv1top                                            6d5s24
              iv1p=iv1+noc(isbv1)                                        6d5s24
              irec=iv1+nvirt(isbv1)*iv2                                 6d5s24
              itri=((iv2*(iv2+isub))/2)+iv1                             6d5s24
              itri=itri+isw*(irec-itri)                                 6d5s24
              sum=0d0                                                   5d31s24
              iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri              6d11s24
              do i=0,naddr(jtv12,ist)-1                                 5d31s24
               sum=sum+bc(jad+i)*bc(iad+i)                                  5d31s24
              end do                                                    5d31s24
              sum=sum*phase                                             6d11s24
              do ii=0,nan(lsa)-1                                        5d31s24
               iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)               5d31s24
               iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv2                     6d11s24
     $                    +nbasisp(jsbv2)*iv1p))                        6d11s24
               iadtt=itt4v(lsa,3)+iv2+nvirt(isbv2)*iq                   6d7s24
               bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                      5d31s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          call second(time2)
          telapo(3)=telapo(3)+time2-time1-tovr
          else                                                          8d26s24
          call second(time2)
          end if                                                        8d26s24
          if(itt4v(lsa,7).gt.0)then                                     8d26s24
          nrow=nbn(lsb)*nbasisp(jsbv2)                                  8d23s24
          itmp=ibcoff                                                   8d23s24
          isum=itmp+nrow*naddr(jtv12,ist)                               8d23s24
          ibcoff=isum+nvv*nrow                                          8d23s24
          call enough('tt4v.sum33',bc,ibc)                              8d23s24
          jtmp=itmp                                                     8d23s24
          do jj=0,nbn(lsb)-1                                            5d31s24
           jv1=ibc(iptbp(lsb)+jj)                                       6d11s24
           do jv2=0,nbasisp(jsbv2)-1                                    6d11s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do i=0,naddr(jtv12,ist)-1                                   8d23s24
             bc(jtmp+i)=bc(jad+i)                                       8d27s24
            end do                                                      8d23s24
            jtmp=jtmp+naddr(jtv12,ist)                                  8d27s24
           end do                                                       8d23s24
          end do                                                        8d23s24
          call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),phase,           8d27s24
     $        bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0, 8d27s24
     $        bc(isum),nvv,'tt4v.sum33')                                8d27s24
          nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv2)                        8d23s24
          do jj=0,nbn(lsb)-1                                            5d31s24
           jv1=ibc(iptbp(lsb)+jj)                                       6d11s24
           do ii=0,nan(lsa)-1                                           8d23s24
            iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)                  8d23s24
            iadtt=itt4v(lsa,7)+nvirt(isbv2)*iq                          8d23s24
            do jv2=0,nbasisp(jsbv2)-1                                    6d11s24
             lsum=isum+nvv*(jv2+nbasisp(jsbv2)*jj)                      8d27s24
             iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv2                       8d23s24
     $                    +nbasisp(jsbv2)*noc(isbv1)))                  8d23s24
             do iv2=0,nvirt(isbv2)-1                                     6d6s24
              iv1top=iv2+1-ist                                           6d5s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
              irec=nvirt(isbv1)*iv2                                     8d23s24
              itri=((iv2*(iv2+isub))/2)                                 8d23s24
              itri=lsum+itri+isw*(irec-itri)                            8d27s24
              do iv1=0,iv1top                                            6d5s24
               bc(iadtt+iv2)=bc(iadtt+iv2)                              8d23s24
     $              +bc(itri+iv1)*eraw2(iadv+nrowi*iv1)                 8d27s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          ibcoff=itmp                                                   8d23s24
          call second(time3)
          telapo(7)=telapo(7)+time3-time2-tovr
          end if                                                        8d26s24
         end if                                                         5d31s24
         isub=-1                                                        5d31s24
         nvv=nvvt                                                       5d31s24
         phase=-1d0                                                     6d11s24
        end do                                                          5d31s24
       end if                                                           5d31s24
      end if                                                            5d31s24
      if(.not.lsame)then                                                5d31s24
       isbv2=lsb                                                        5d31s24
       isbv1=lsd                                                        5d31s24
       jsbv1=lsa                                                        6d11s24
       jsbv2=lsc                                                        6d11s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         5d31s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
        if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
         nv1=nvirt(isbv1)                                               6d5s24
         nv2=nvirt(isbv2)                                               6d5s24
         mv1=nbasisp(jsbv1)                                             6d10s24
         mv2=nbasisp(jsbv2)                                             6d10s24
         if(isbv2.eq.isbv1)then                                         6d7s24
          isw=0                                                          5d20s24
          nvvs=(nv1*(nv1+1))/2                                          5d31s24
          nvvt=(nv1*(nv1-1))/2                                          5d31s24
         else                                                            5d20s24
          isw=1                                                          5d20s24
          nvvs=nv1*nv2                                                  5d31s24
          nvvt=nvvs                                                      5d20s24
         end if                                                          5d20s24
         isub=1                                                          5d31s24
         nvv=nvvs                                                        5d31s24
         mvv=mv1*mv2                                                    6d7s24
         phase=1d0                                                      6d11s24
         do ist=1,2                                                      5d31s24
          if(min(nvv,naddr(jtv12,ist)).gt.0)then                        8d23s24
           if(itt4v(lsb,3).gt.0)then                                    8d26s24
           call second(time1)
           do jj=0,nan(lsa)-1                                            5d31s24
            jv1=ibc(iptap(lsa)+jj)                                      6d11s24
            do jv2=0,nbasisp(jsbv2)-1                                   6d11s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do iv2=0,nvirt(isbv2)-1                                    6d5s24
              iv1top=iv2+1-ist                                           5d31s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
              do iv1=0,iv1top                                            5d31s24
               iv1p=iv1+noc(isbv1)                                       6d5s24
               irec=iv1+nvirt(isbv1)*iv2                                6d5s24
               itri=((iv2*(iv2+isub))/2)+iv1                             5d31s24
               itri=itri+isw*(irec-itri)                                 5d31s24
               sum=0d0                                                  5d31s24
               iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri             6d11s24
               do i=0,naddr(jtv12,ist)-1                                5d31s24
                sum=sum+bc(jad+i)*bc(iad+i)                             6d11s24
               end do                                                   5d31s24
               sum=sum*phase                                            6d11s24
               do ii=0,nbn(lsb)-1                                       5d31s24
                iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)              5d31s24
                iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv2                    6d11s24
     $                    +nbasisp(jsbv2)*iv1p))                        6d11s24
                iadtt=itt4v(lsb,3)+iv2+nvirt(isbv2)*iq                  6d7s24
                bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                     5d31s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           call second(time2)
           telapo(3)=telapo(3)+time2-time1-tovr
           else                                                         8d26s24
           call second(time2)
           end if                                                       8d26s24
           if(itt4v(lsb,7).gt.0)then                                    8d26s24
           nrow=nan(lsa)*nbasisp(jsbv2)                                 8d23s24
           itmp=ibcoff                                                  8d23s24
           isum=itmp+nrow*naddr(jtv12,ist)                              8d23s24
           ibcoff=isum+nrow*nvv                                         8d23s24
           call enough('tt4v.sum34',bc,ibc)                             8d23s24
           jtmp=itmp                                                    8d23s24
           do jj=0,nan(lsa)-1                                            5d31s24
            jv1=ibc(iptap(lsa)+jj)                                      6d11s24
            do jv2=0,nbasisp(jsbv2)-1                                   6d11s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do i=0,naddr(jtv12,ist)-1                                  8d23s24
              bc(jtmp+i)=bc(jad+i)                                      8d27s24
             end do                                                     8d23s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d27s24
            end do                                                      8d23s24
           end do                                                       8d23s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),phase,          8d27s24
     $       bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,  8d27s24
     $       bc(isum),nvv,'tt4v.sum34')                                 8d27s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv2)                       8d23s24
           do jj=0,nan(lsa)-1                                            5d31s24
            jv1=ibc(iptap(lsa)+jj)                                      6d11s24
            do ii=0,nbn(lsb)-1                                          8d23s24
             iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)                 8d23s24
             iadtt=itt4v(lsb,7)+nvirt(isbv2)*iq                         8d23s24
             do jv2=0,nbasisp(jsbv2)-1                                   6d11s24
              lsum=isum+nvv*(jv2+nbasisp(jsbv2)*jj)                     8d27s24
              iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv2                      8d23s24
     $                    +nbasisp(jsbv2)*noc(isbv1)))                  8d23s24
              do iv2=0,nvirt(isbv2)-1                                    6d5s24
               iv1top=iv2+1-ist                                           5d31s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
               irec=nvirt(isbv1)*iv2                                    8d23s24
               itri=((iv2*(iv2+isub))/2)                                8d23s24
               itri=lsum+itri+isw*(irec-itri)                           8d27s24
               do iv1=0,iv1top                                            5d31s24
                bc(iadtt+iv2)=bc(iadtt+iv2)                             8d23s24
     $               +bc(itri+iv1)*eraw2(iadv+nrowi*iv1)                8d27s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           ibcoff=itmp                                                  8d23s24
           call second(time3)
           telapo(7)=telapo(7)+time3-time2-tovr
           end if                                                       8d26s24
          end if                                                         5d31s24
          isub=-1                                                       5d31s24
          nvv=nvvt                                                      5d31s24
          phase=-1d0                                                    6d11s24
         end do                                                          5d31s24
        end if                                                           5d31s24
       end if                                                            5d31s24
      end if                                                            5d31s24
cc    k4
      isbv2=lsb                                                         6d10s24
      isbv1=lsd                                                         6d10s24
      jsbv1=lsa                                                         6d11s24
      jsbv2=lsc                                                         6d11s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d6s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
       if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
        nv1=nvirt(isbv1)                                                6d5s24
        nv2=nvirt(isbv2)                                                6d5s24
        mv1=nbasisp(jsbv1)                                              6d10s24
        mv2=nbasisp(jsbv2)                                              6d10s24
        if(isbv2.eq.isbv1)then                                          6d7s24
         isw=0                                                          5d20s24
         nvvs=(nv1*(nv1+1))/2                                           5d31s24
         nvvt=(nv1*(nv1-1))/2                                           5d31s24
        else                                                            5d20s24
         isw=1                                                          5d20s24
         nvvs=nv1*nv2                                                   5d31s24
         nvvt=nvvs                                                      5d20s24
        end if                                                          5d20s24
        isub=1                                                          5d31s24
        nvv=nvvs                                                        5d31s24
        mvv=mv1*mv2                                                     6d7s24
        phase=1d0                                                       6d11s24
        do ist=1,2                                                      5d31s24
         if(min(nvv,naddr(jtv12,ist)).gt.0)then                         8d23s24
          if(itt4v(lsb,4).gt.0)then                                     8d26s24
          call second(time1)
          do jj=0,nan(lsa)-1                                            5d31s24
           jv1=ibc(iptap(lsa)+jj)                                       6d11s24
           do jv2=0,nbasisp(jsbv2)-1                                    6d11s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do iv2=0,nvirt(isbv2)-1                                     6d6s24
             iv1top=iv2+1-ist                                           6d5s24
             iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
             do iv1=0,iv1top                                            6d5s24
              iv1p=iv1+noc(isbv1)                                       6d10s24
              irec=iv1+nvirt(isbv1)*iv2                                 6d5s24
              itri=((iv2*(iv2+isub))/2)+iv1                             6d5s24
              itri=itri+isw*(irec-itri)                                 6d5s24
              sum=0d0                                                   5d31s24
              iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri              6d11s24
              do i=0,naddr(jtv12,ist)-1                                 5d31s24
               sum=sum+bc(jad+i)*bc(iad+i)                              6d11s24
              end do                                                    5d31s24
              sum=sum*phase                                             6d11s24
              do ii=0,nbn(lsb)-1                                        5d31s24
               iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)               5d31s24
               iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv2                     6d11s24
     $                    +nbasisp(jsbv2)*iv1p))                        6d11s24
               iadtt=itt4v(lsb,4)+iv2+nvirt(isbv2)*iq                   6d10s24
               bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                      5d31s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          call second(time2)
          telapo(4)=telapo(4)+time2-time1-tovr
          else                                                          8d26s24
          call second(time2)
          end if                                                        8d26s24
          if(itt4v(lsb,8).gt.0)then                                     8d26s24
          nrow=nan(lsa)*nbasisp(jsbv2)                                  8d23s24
          itmp=ibcoff                                                   8d23s24
          isum=itmp+nrow*naddr(jtv12,ist)                               8d23s24
          ibcoff=isum+nrow*nvv                                          8d23s24
          call enough('tt4v.sum43',bc,ibc)                              8d23s24
          jtmp=itmp                                                     8d23s24
          do jj=0,nan(lsa)-1                                            5d31s24
           jv1=ibc(iptap(lsa)+jj)                                       6d11s24
           do jv2=0,nbasisp(jsbv2)-1                                    6d11s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
            jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec             6d11s24
            do i=0,naddr(jtv12,ist)-1                                   8d23s24
             bc(jtmp+i)=bc(jad+i)                                       8d23s24
            end do                                                      8d23s24
            jtmp=jtmp+naddr(jtv12,ist)                                  8d23s24
           end do                                                       8d23s24
          end do                                                        8d23s24
          call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),phase,           8d23s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d23s24
     $         bc(isum),nvv,'tt4v.sum43')                               8d23s24
          nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv2)                        8d23s24
          do jj=0,nan(lsa)-1                                            5d31s24
           jv1=ibc(iptap(lsa)+jj)                                       6d11s24
           do ii=0,nbn(lsb)-1                                           8d23s24
            iq=ibc(iptbp(lsb)+ii)+iqshift*nbasisp(lsb)                  8d23s24
            iadtt=itt4v(lsb,8)+nvirt(isbv2)*iq                          8d23s24
            do jv2=0,nbasisp(jsbv2)-1                                    6d11s24
             lsum=isum+nvv*(jv2+nbasisp(jsbv2)*jj)                      8d23s24
             iadv=1+jj+nan(lsa)*(ii+nbn(lsb)*(jv2                       8d23s24
     $                    +nbasisp(jsbv2)*noc(isbv1)))                  8d23s24
             do iv2=0,nvirt(isbv2)-1                                     6d6s24
              iv1top=iv2+1-ist                                           6d5s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                  6d6s24
              irec=nvirt(isbv1)*iv2                                     8d23s24
              itri=((iv2*(iv2+isub))/2)                                 8d23s24
              itri=lsum+itri+isw*(irec-itri)                            8d23s24
              do iv1=0,iv1top                                            6d5s24
               bc(iadtt+iv2)=bc(iadtt+iv2)+bc(itri+iv1)                 8d23s24
     $              *eraw2(iadv+nrowi*iv1)                              8d23s24
              end do                                                    5d31s24
             end do                                                     5d31s24
            end do                                                      5d31s24
           end do                                                       5d31s24
          end do                                                        5d31s24
          ibcoff=itmp                                                   8d23s24
          call second(time3)
          telapo(8)=telapo(8)+time3-time2-tovr
          end if                                                        8d26s24
         end if                                                         5d31s24
         isub=-1                                                        5d31s24
         nvv=nvvt                                                       5d31s24
         phase=-1d0                                                     6d11s24
        end do                                                          5d31s24
       end if                                                           5d31s24
      end if                                                            5d31s24
      if(.not.lsame)then                                                5d31s24
       isbv2=lsa                                                        6d10s24
       isbv1=lsd                                                        6d10s24
       jsbv1=lsb                                                        6d11s24
       jsbv2=lsc                                                        6d11s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         5d31s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                 5d31s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                 5d31s24
        if(naddr(itv12,3).eq.naddr(jtv12,3).and.                           5d31s24
     $      max(naddr(itv12,1),naddr(itv12,2)).gt.0)then                  5d31s24
         nv1=nvirt(isbv1)                                               6d5s24
         nv2=nvirt(isbv2)                                               6d5s24
         mv1=nbasisp(jsbv1)                                             6d10s24
         mv2=nbasisp(jsbv2)                                             6d10s24
         if(isbv2.eq.isbv1)then                                         6d7s24
          isw=0                                                          5d20s24
          nvvs=(nv1*(nv1+1))/2                                          5d31s24
          nvvt=(nv1*(nv1-1))/2                                          5d31s24
         else                                                            5d20s24
          isw=1                                                          5d20s24
          nvvs=nv1*nv2                                                  5d31s24
          nvvt=nvvs                                                      5d20s24
         end if                                                          5d20s24
         isub=1                                                          5d31s24
         nvv=nvvs                                                        5d31s24
         mvv=mv1*mv2                                                    6d7s24
         phase=1d0                                                      6d11s24
         do ist=1,2                                                      5d31s24
          if(min(nvv,naddr(jtv12,ist)).gt.0)then                        8d23s24
           if(itt4v(lsa,4).gt.0)then                                    8d26s24
           call second(time1)
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv1=ibc(iptbp(lsb)+jj)                                      6d11s24
            do jv2=0,nbasisp(jsbv2)-1                                   6d11s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do iv2=0,nvirt(isbv2)-1                                    6d5s24
              iv1top=iv2+1-ist                                           5d31s24
              iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
              do iv1=0,iv1top                                            5d31s24
               iv1p=iv1+noc(isbv1)                                      6d10s24
               irec=iv1+nvirt(isbv1)*iv2                                6d5s24
               itri=((iv2*(iv2+isub))/2)+iv1                             5d31s24
               itri=itri+isw*(irec-itri)                                 5d31s24
               sum=0d0                                                  5d31s24
               iad=iaddr(itv12,ist,5)+naddr(jtv12,ist)*itri             6d11s24
               do i=0,naddr(jtv12,ist)-1                                5d31s24
                sum=sum+bc(jad+i)*bc(iad+i)                             6d11s24
               end do                                                   5d31s24
               sum=sum*phase                                            6d11s24
               do ii=0,nan(lsa)-1                                       5d31s24
                iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)              5d31s24
                iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv2                    6d11s24
     $                    +nbasisp(jsbv2)*iv1p))                        6d11s24
                iadtt=itt4v(lsa,4)+iv2+nvirt(isbv2)*iq                  6d10s24
                bc(iadtt)=bc(iadtt)+sum*eraw2(iadv)                     5d31s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           call second(time2)
           telapo(4)=telapo(4)+time2-time1-tovr
           else                                                         8d26s24
           call second(time2)
           end if                                                       8d26s24
           if(itt4v(lsa,8).gt.0)then                                    8d26s24
           nrow=nbn(lsb)*nbasisp(jsbv2)                                 8d23s24
           itmp=ibcoff                                                  8d23s24
           isum=itmp+nrow*naddr(jtv12,ist)                              8d23s24
           ibcoff=isum+nrow*nvv                                         8d23s24
           call enough('tt4b.sum44',bc,ibc)                             8d23s24
           jtmp=itmp                                                    8d23s24
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv1=ibc(iptbp(lsb)+jj)                                      6d11s24
            do jv2=0,nbasisp(jsbv2)-1                                   6d11s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 5d31s24
             jad=iaddr(jtv12,ist,ier2)+naddr(jtv12,ist)*jrec            6d11s24
             do i=0,naddr(jtv12,ist)-1                                  8d23s24
              bc(jtmp+i)=bc(jad+i)                                      8d23s24
             end do                                                     8d23s24
             jtmp=jtmp+naddr(jtv12,ist)                                 8d23s24
            end do                                                      8d23s24
           end do                                                       8d23s24
           call dgemm('n','n',nvv,nrow,naddr(jtv12,ist),phase,          8d23s24
     $         bc(iaddr(itv12,ist,6)),nvv,bc(itmp),naddr(jtv12,ist),0d0,8d23s24
     $          bc(isum),nvv,'tt4v.sum44')                              8d23s24
           nrowi=nan(lsa)*nbn(lsb)*nbasisp(jsbv2)                       8d23s24
           do jj=0,nbn(lsb)-1                                            5d31s24
            jv1=ibc(iptbp(lsb)+jj)                                      6d11s24
            do ii=0,nan(lsa)-1                                          8d23s24
             iq=ibc(iptap(lsa)+ii)+iqshift*nbasisp(lsa)                 8d23s24
             iadtt=itt4v(lsa,8)+nvirt(isbv2)*iq                         8d23s24
             do jv2=0,nbasisp(jsbv2)-1                                   6d11s24
              lsum=isum+nvv*(jv2+nbasisp(jsbv2)*jj)                     8d23s24
              iadv=1+ii+nan(lsa)*(jj+nbn(lsb)*(jv2                      8d23s24
     $                    +nbasisp(jsbv2)*noc(isbv1)))                  8d23s24
              do iv2=0,nvirt(isbv2)-1                                    6d5s24
               iv1top=iv2+1-ist                                           5d31s24
               iv1top=iv1top+isw*(nvirt(isbv1)-1-iv1top)                 6d5s24
               irec=nvirt(isbv1)*iv2                                    8d23s24
               itri=((iv2*(iv2+isub))/2)                                8d23s24
               itri=lsum+itri+isw*(irec-itri)                           8d23s24
               do iv1=0,iv1top                                            5d31s24
                bc(iadtt+iv2)=bc(iadtt+iv2)                             8d23s24
     $               +bc(itri+iv1)*eraw2(iadv+nrowi*iv1)                8d23s24
               end do                                                   5d31s24
              end do                                                     5d31s24
             end do                                                      5d31s24
            end do                                                       5d31s24
           end do                                                        5d31s24
           ibcoff=itmp                                                  8d23s24
           call second(time3)
           telapo(8)=telapo(8)+time3-time2-tovr
           end if                                                       8d26s24
          end if                                                         5d31s24
          isub=-1                                                       5d31s24
          nvv=nvvt                                                      5d31s24
          phase=-1d0                                                    6d11s24
         end do                                                          5d31s24
        end if                                                           5d31s24
       end if                                                            5d31s24
      end if                                                            5d31s24
      return                                                            5d20s24
      end                                                               5d20s24
