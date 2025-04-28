c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine der4v(iaddr,naddr,d4v,lsa,lsb,lsc,lsd,nsymb,multh,     6d13s24
     $     nbasisp,nerimat,nan,iptap,nbn,iptbp,in,ncomp,lsame,derint,   6d13s24
     $     bc,ibc)                                                      6d13s24
      implicit real*8 (a-h,o-z)                                         6d13s24
      include "common.store"                                            6d13s24
      logical lsame                                                     6d13s24
      dimension iaddr(36,2,4),naddr(36,3),multh(8,8),nbasisp(*),nan(*), 6d13s24
     $     iptap(*),nbn(*),iptbp(*),derint(*),d4v(4)                    6d13s24
      data icall/0/
      save
      icall=icall+1
      nxnx=nan(lsa)*nbn(lsb)*nbasisp(lsc)*nbasisp(lsd)
      ifound=0
c
c     in           vec
c     0  LLLL      LL
c     1  LLSS      LS
c     2  SSLL      SL
c     3  SSSS      SS
c
      ier=1                                                             6d13s24
      jer=1                                                             6d13s24
      isbv1=lsa                                                         6d13s24
      jsbv1=lsb                                                         6d13s24
      isbv2=lsc                                                         6d13s24
      jsbv2=lsd                                                         6d13s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d13s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                6d13s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                6d13s24
       if(naddr(itv12,3).eq.naddr(jtv12,3))then                         6d13s24
        if(in.eq.0)then                                                 6d14s24
         ier=1                                                          6d14s24
         jer=1                                                          6d14s24
        else if(in.eq.1)then                                            6d14s24
         ier=2                                                          6d14s24
         jer=2                                                          6d14s24
        else if(in.eq.2)then                                            6d14s24
         ier=3                                                          6d14s24
         jer=3                                                          6d14s24
        else                                                            6d14s24
         ier=4                                                          6d14s24
         jer=4                                                          6d14s24
        end if                                                          6d14s24
        if(isbv1.eq.isbv2.and.(in.eq.0.or.in.eq.3))then                 6d13s24
         isw=0                                                          6d13s24
         nvvs=(nbasisp(isbv1)*(nbasisp(isbv1)+1))/2                     6d13s24
         mvvs=(nbasisp(jsbv1)*(nbasisp(jsbv1)+1))/2                     6d13s24
         nvvt=(nbasisp(isbv1)*(nbasisp(isbv1)-1))/2                     6d13s24
         mvvt=(nbasisp(jsbv1)*(nbasisp(jsbv1)-1))/2                     6d13s24
        else                                                            6d13s24
         isw=1                                                          6d13s24
         nvvs=nbasisp(isbv1)*nbasisp(isbv2)                             6d13s24
         nvvt=nvvs                                                      6d13s24
         mvvs=nbasisp(jsbv1)*nbasisp(jsbv2)                             6d13s24
         mvvt=mvvs                                                      6d13s24
        end if                                                          6d13s24
        nvv=nvvs                                                        6d13s24
        mvv=mvvs                                                        6d13s24
        isub=+1                                                         6d13s24
        do ist=1,2                                                      6d13s24
         do jj=0,nbn(lsb)-1                                             6d13s24
          jv1=ibc(iptbp(lsb)+jj)                                        6d13s24
          jv2bot=jv1+ist-1                                              6d13s24
          jv2bot=jv2bot*(1-isw)                                         6d13s24
          do ii=0,nan(lsa)-1                                            6d13s24
           iv1=ibc(iptap(lsa)+ii)                                       6d13s24
           iv2bot=iv1+ist-1                                             6d13s24
           iv2bot=iv2bot*(1-isw)                                        6d13s24
           do jv2=jv2bot,nbasisp(jsbv2)-1                                6d13s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 6d13s24
            jtri=((jv2*(jv2+isub))/2)+jv1                               6d13s24
            jtri=jtri+isw*(jrec-jtri)                                   6d13s24
            do iv2=iv2bot,nbasisp(isbv2)-1                              6d13s24
             irec=iv1+nbasisp(isbv1)*iv2                                6d13s24
             itri=((iv2*(iv2+isub))/2)+iv1                              6d13s24
             itri=itri+isw*(irec-itri)                                  6d13s24
             sum=0d0                                                    6d13s24
             iad=iaddr(itv12,ist,ier)+naddr(itv12,ist)*itri             6d13s24
             jad=iaddr(jtv12,ist,jer)+naddr(itv12,ist)*jtri             6d13s24
             do i=0,naddr(itv12,ist)-1                                  6d13s24
              sum=sum+bc(iad+i)*bc(jad+i)                               6d13s24
             end do                                                     6d13s24
             iadd=ii+1+nan(lsa)*(jj+nbn(lsb)*(iv2+nbasisp(lsc)*jv2))    6d13s24
             orig=d4v(1)
             d4v(1)=d4v(1)+sum*derint(iadd)                             6d13s24
            end do                                                      6d13s24
           end do                                                       6d13s24
          end do                                                          6d13s24
         end do                                                         6d13s24
         nvv=nvvt                                                       6d13s24
         mvv=mvvt                                                       6d13s24
         isub=-1                                                        6d13s24
        end do                                                          6d13s24
       end if                                                           6d13s24
      end if                                                            6d13s24
      isbv1=lsa                                                         6d13s24
      jsbv2=lsb                                                         6d13s24
      isbv2=lsc                                                         6d13s24
      jsbv1=lsd                                                         6d13s24
      if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d13s24
       itv12=((isbv2*(isbv2-1))/2)+isbv1                                6d13s24
       jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                6d13s24
       if(naddr(itv12,3).eq.naddr(jtv12,3))then                         6d13s24
        if(in.eq.0)then                                                 6d14s24
         ier=1                                                          6d14s24
         jer=1                                                          6d14s24
        else if(in.eq.1)then                                            6d14s24
         ier=2                                                          6d14s24
         jer=3                                                          6d14s24
        else if(in.eq.2)then                                            6d14s24
         ier=3                                                          6d14s24
         jer=2                                                          6d14s24
        else                                                            6d14s24
         ier=4                                                          6d14s24
         jer=4                                                          6d14s24
        end if                                                          6d14s24
        if(isbv1.eq.isbv2.and.(in.eq.0.or.in.eq.3))then                 6d13s24
         isw=0                                                          6d13s24
         nvvs=(nbasisp(isbv1)*(nbasisp(isbv1)+1))/2                     6d13s24
         mvvs=(nbasisp(jsbv1)*(nbasisp(jsbv1)+1))/2                     6d13s24
         nvvt=(nbasisp(isbv1)*(nbasisp(isbv1)-1))/2                     6d13s24
         mvvt=(nbasisp(jsbv1)*(nbasisp(jsbv1)-1))/2                     6d13s24
        else                                                            6d13s24
         isw=1                                                          6d13s24
         nvvs=nbasisp(isbv1)*nbasisp(isbv2)                             6d13s24
         nvvt=nvvs                                                      6d13s24
         mvvs=nbasisp(jsbv1)*nbasisp(jsbv2)                             6d13s24
         mvvt=mvvs                                                      6d13s24
        end if                                                          6d13s24
        nvv=nvvs                                                        6d13s24
        mvv=mvvs                                                        6d13s24
        isub=+1                                                         6d13s24
        phase=1d0                                                       6d13s24
        do ist=1,2                                                      6d13s24
         do jj=0,nbn(lsb)-1                                             6d13s24
          jv2=ibc(iptbp(lsb)+jj)                                        6d13s24
          jv1top=jv2+1-ist                                              6d13s24
          jv1top=jv1top+isw*(nbasisp(jsbv1)-1-jv1top)                   6d13s24
          do ii=0,nan(lsa)-1                                            6d13s24
           iv1=ibc(iptap(lsa)+ii)                                       6d13s24
           iv2bot=iv1+ist-1                                             6d13s24
           iv2bot=iv2bot*(1-isw)                                        6d13s24
           do jv1=0,jv1top                                              6d13s24
            jrec=jv1+nbasisp(jsbv1)*jv2                                 6d13s24
            jtri=((jv2*(jv2+isub))/2)+jv1                               6d13s24
            jtri=jtri+isw*(jrec-jtri)                                   6d13s24
            do iv2=iv2bot,nbasisp(isbv2)-1                              6d13s24
             irec=iv1+nbasisp(isbv1)*iv2                                6d13s24
             itri=((iv2*(iv2+isub))/2)+iv1                              6d13s24
             itri=itri+isw*(irec-itri)                                  6d13s24
             sum=0d0                                                    6d13s24
             iad=iaddr(itv12,ist,ier)+naddr(itv12,ist)*itri             6d13s24
             jad=iaddr(jtv12,ist,jer)+naddr(itv12,ist)*jtri             6d13s24
             do i=0,naddr(itv12,ist)-1                                  6d13s24
              sum=sum+bc(iad+i)*bc(jad+i)                               6d13s24
             end do                                                     6d13s24
             sum=sum*phase                                              6d13s24
             iadd=ii+1+nan(lsa)*(jj+nbn(lsb)*(iv2+nbasisp(lsc)*jv1))    6d13s24
             orig=d4v(2)
             d4v(2)=d4v(2)+sum*derint(iadd)                             6d13s24
            end do                                                      6d13s24
           end do                                                       6d13s24
          end do                                                          6d13s24
         end do                                                         6d13s24
         nvv=nvvt                                                       6d13s24
         mvv=mvvt                                                       6d13s24
         isub=-1                                                        6d13s24
         phase=-1d0                                                     6d13s24
        end do                                                          6d13s24
       end if                                                           6d13s24
      end if                                                            6d13s24
      if(.not.lsame)then                                                6d13s24
       isbv1=lsb                                                         6d13s24
       jsbv1=lsa                                                         6d13s24
       isbv2=lsc                                                         6d13s24
       jsbv2=lsd                                                         6d13s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d13s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                6d13s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                6d13s24
        if(naddr(itv12,3).eq.naddr(jtv12,3))then                         6d13s24
         if(in.eq.0)then                                                 6d14s24
          ier=1                                                          6d14s24
          jer=1                                                          6d14s24
         else if(in.eq.1)then                                            6d14s24
          ier=2                                                          6d14s24
          jer=2                                                          6d14s24
         else if(in.eq.2)then                                            6d14s24
          ier=3                                                          6d14s24
          jer=3                                                          6d14s24
         else                                                            6d14s24
          ier=4                                                          6d14s24
          jer=4                                                          6d14s24
         end if                                                          6d14s24
         if(isbv1.eq.isbv2.and.(in.eq.0.or.in.eq.3))then                 6d13s24
          isw=0                                                          6d13s24
          nvvs=(nbasisp(isbv1)*(nbasisp(isbv1)+1))/2                     6d13s24
          mvvs=(nbasisp(jsbv1)*(nbasisp(jsbv1)+1))/2                     6d13s24
          nvvt=(nbasisp(isbv1)*(nbasisp(isbv1)-1))/2                     6d13s24
          mvvt=(nbasisp(jsbv1)*(nbasisp(jsbv1)-1))/2                     6d13s24
         else                                                            6d13s24
          isw=1                                                          6d13s24
          nvvs=nbasisp(isbv1)*nbasisp(isbv2)                             6d13s24
          nvvt=nvvs                                                      6d13s24
          mvvs=nbasisp(jsbv1)*nbasisp(jsbv2)                             6d13s24
          mvvt=mvvs                                                      6d13s24
         end if                                                          6d13s24
         nvv=nvvs                                                        6d13s24
         mvv=mvvs                                                        6d13s24
         isub=+1                                                         6d13s24
         do ist=1,2                                                      6d13s24
          do jj=0,nbn(lsb)-1                                            6d13s24
           iv1=ibc(iptbp(lsb)+jj)                                       6d13s24
           iv2bot=iv1+ist-1                                             6d13s24
           iv2bot=iv2bot*(1-isw)                                        6d13s24
           do ii=0,nan(lsa)-1                                           6d13s24
            jv1=ibc(iptap(lsa)+ii)                                      6d13s24
            jv2bot=jv1+ist-1                                            6d13s24
            jv2bot=jv2bot*(1-isw)                                       6d13s24
            do jv2=jv2bot,nbasisp(jsbv2)-1                                6d13s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 6d13s24
             jtri=((jv2*(jv2+isub))/2)+jv1                               6d13s24
             jtri=jtri+isw*(jrec-jtri)                                   6d13s24
             do iv2=iv2bot,nbasisp(isbv2)-1                              6d13s24
              irec=iv1+nbasisp(isbv1)*iv2                                6d13s24
              itri=((iv2*(iv2+isub))/2)+iv1                              6d13s24
              itri=itri+isw*(irec-itri)                                  6d13s24
              sum=0d0                                                    6d13s24
              iad=iaddr(itv12,ist,ier)+naddr(itv12,ist)*itri             6d13s24
              jad=iaddr(jtv12,ist,jer)+naddr(itv12,ist)*jtri             6d13s24
              do i=0,naddr(itv12,ist)-1                                  6d13s24
               sum=sum+bc(iad+i)*bc(jad+i)                               6d13s24
              end do                                                     6d13s24
              iadd=ii+1+nan(lsa)*(jj+nbn(lsb)*(iv2+nbasisp(lsc)*jv2))    6d13s24
              orig=d4v(3)
              d4v(3)=d4v(3)+sum*derint(iadd)                            6d13s24
             end do                                                      6d13s24
            end do                                                       6d13s24
           end do                                                          6d13s24
          end do                                                         6d13s24
          nvv=nvvt                                                       6d13s24
          mvv=mvvt                                                       6d13s24
          isub=-1                                                        6d13s24
         end do                                                          6d13s24
        end if                                                           6d13s24
       end if                                                            6d13s24
       isbv1=lsb                                                         6d13s24
       jsbv2=lsa                                                         6d13s24
       isbv2=lsc                                                         6d13s24
       jsbv1=lsd                                                         6d13s24
       if(isbv2.ge.isbv1.and.jsbv2.ge.jsbv1)then                         6d13s24
        itv12=((isbv2*(isbv2-1))/2)+isbv1                                6d13s24
        jtv12=((jsbv2*(jsbv2-1))/2)+jsbv1                                6d13s24
        if(naddr(itv12,3).eq.naddr(jtv12,3))then                         6d13s24
         if(in.eq.0)then                                                 6d14s24
          ier=1                                                          6d14s24
          jer=1                                                          6d14s24
         else if(in.eq.1)then                                            6d14s24
          ier=2                                                          6d14s24
          jer=3                                                          6d14s24
         else if(in.eq.2)then                                            6d14s24
          ier=3                                                          6d14s24
          jer=2                                                          6d14s24
         else                                                            6d14s24
          ier=4                                                          6d14s24
          jer=4                                                          6d14s24
         end if                                                          6d14s24
         if(isbv1.eq.isbv2.and.(in.eq.0.or.in.eq.3))then                 6d13s24
          isw=0                                                          6d13s24
          nvvs=(nbasisp(isbv1)*(nbasisp(isbv1)+1))/2                     6d13s24
          mvvs=(nbasisp(jsbv1)*(nbasisp(jsbv1)+1))/2                     6d13s24
          nvvt=(nbasisp(isbv1)*(nbasisp(isbv1)-1))/2                     6d13s24
          mvvt=(nbasisp(jsbv1)*(nbasisp(jsbv1)-1))/2                     6d13s24
         else                                                            6d13s24
          isw=1                                                          6d13s24
          nvvs=nbasisp(isbv1)*nbasisp(isbv2)                             6d13s24
          nvvt=nvvs                                                      6d13s24
          mvvs=nbasisp(jsbv1)*nbasisp(jsbv2)                             6d13s24
          mvvt=mvvs                                                      6d13s24
         end if                                                          6d13s24
         nvv=nvvs                                                        6d13s24
         mvv=mvvs                                                        6d13s24
         isub=+1                                                         6d13s24
         phase=1d0                                                       6d13s24
         do ist=1,2                                                      6d13s24
          do jj=0,nbn(lsb)-1                                             6d13s24
           iv1=ibc(iptbp(lsb)+jj)                                        6d13s24
           iv2bot=iv1+ist-1                                             6d13s24
           iv2bot=iv2bot*(1-isw)                                        6d13s24
           do ii=0,nan(lsa)-1                                           6d13s24
            jv2=ibc(iptap(lsa)+ii)                                      6d13s24
            jv1top=jv2+1-ist                                              6d13s24
            jv1top=jv1top+isw*(nbasisp(jsbv1)-1-jv1top)                   6d13s24
            do jv1=0,jv1top                                              6d13s24
             jrec=jv1+nbasisp(jsbv1)*jv2                                 6d13s24
             jtri=((jv2*(jv2+isub))/2)+jv1                               6d13s24
             jtri=jtri+isw*(jrec-jtri)                                   6d13s24
             do iv2=iv2bot,nbasisp(isbv2)-1                              6d13s24
              irec=iv1+nbasisp(isbv1)*iv2                                6d13s24
              itri=((iv2*(iv2+isub))/2)+iv1                              6d13s24
              itri=itri+isw*(irec-itri)                                  6d13s24
              sum=0d0                                                    6d13s24
              iad=iaddr(itv12,ist,ier)+naddr(itv12,ist)*itri             6d13s24
              jad=iaddr(jtv12,ist,jer)+naddr(itv12,ist)*jtri             6d13s24
              do i=0,naddr(itv12,ist)-1                                  6d13s24
               sum=sum+bc(iad+i)*bc(jad+i)                               6d13s24
              end do                                                     6d13s24
              sum=sum*phase                                              6d13s24
              iadd=ii+1+nan(lsa)*(jj+nbn(lsb)*(iv2+nbasisp(lsc)*jv1))    6d13s24
              orig=d4v(4)
              d4v(4)=d4v(4)+sum*derint(iadd)                            6d13s24
             end do                                                      6d13s24
            end do                                                       6d13s24
           end do                                                          6d13s24
          end do                                                         6d13s24
          nvv=nvvt                                                       6d13s24
          mvv=mvvt                                                       6d13s24
          isub=-1                                                        6d13s24
          phase=-1d0                                                     6d13s24
         end do                                                          6d13s24
        end if                                                           6d13s24
       end if                                                            6d13s24
      end if                                                            6d13s24
      return                                                            6d13s24
      end                                                               6d13s24
