c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine fillrup(iwaveb,iwavek,itmom,jjlowb,jjlowk,tm,ndim,     6d4s21
     $     joffb,joffk,ioffb,ioffk,ndimb,tmb,bc,ibc)                    11d14s22
      implicit real*8 (a-h,o-z)                                         6d4s21
      character*6 lab
      integer*1 ipack1(4)                                               6d4s21
      integer*8 itmom(*),joffb(*),ioffb(*),joffk(*),ioffk(*)            6d4s21
      equivalence (ipack1,npack4)                                       6d4s21
      dimension iwaveb(*),iwavek(*),tm(ndim,ndimb,2),tmb(ndimb,*)       1d13s23
      include "common.store"                                            6d4s21
      isb=iwaveb(1)
      nrootb=iwaveb(3)                                                  6d4s21
      npack4=iwaveb(6)                                                  6d4s21
      do l=1,6                                                          6d4s21
       if(iwaveb(13+l).ne.0)then                                        6d4s21
        lab(l:l)=char(iwaveb(13+l))                                     6d4s21
       else                                                             6d4s21
        lab(l:l)=' '                                                    6d4s21
       end if                                                           6d4s21
      end do                                                            6d4s21
      isb=isb-1                                                         6d4s21
      lb=ipack1(3)                                                      6d4s21
      lb2=lb*2
      isk=iwavek(1)
      nrootk=iwavek(3)                                                  6d4s21
      npack4=iwavek(6)                                                  6d4s21
      do l=1,6                                                          6d4s21
       if(iwavek(13+l).ne.0)then                                        6d4s21
        lab(l:l)=char(iwavek(13+l))                                     6d4s21
       else                                                             6d4s21
        lab(l:l)=' '                                                    6d4s21
       end if                                                           6d4s21
      end do                                                            6d4s21
      isk=isk-1                                                         6d4s21
      if(isb.ne.isk)return                                              6d4s21
      lk=ipack1(3)                                                      6d4s21
      lk2=lk*2
      jlowb=iabs(lb2-isb)                                               6d4s21
      jhib=lb2+isb                                                      6d4s21
      jlowk=iabs(lk2-isk)                                               6d4s21
      jhik=lk2+isk                                                      6d4s21
      do ipass=1,3
       if(ipass.eq.3)then                                               6d4s21
        iq2=4                                                           6d4s21
       else                                                             6d4s21
        iq2=2                                                           6d4s21
       end if                                                           6d4s21
       if(itmom(ipass).gt.0)then                                           6d4s21
        if(ipass.eq.1)then                                               6d4s21
         isto=1                                                         6d4s21
        else if(ipass.eq.2)then                                               6d4s21
         isto=1                                                         6d4s21
        else                                                             6d4s21
         isto=2                                                         6d4s21
        end if                                                           6d4s21
        do jk=jlowk,jhik,2                                               6d4s21
         jjk=1+((jk-jjlowk)/2)                                            6d4s21
         iket=joffk(jjk)+ioffk(jjk)+1                                   6d4s21
         do jb=jlowb,jhib,2                                              6d4s21
          jjb=1+((jb-jjlowb)/2)                                           6d4s21
          fact=f6j(jb,jk,iq2,lk2,lb2,isb,1)                              6d4s21
          fact=fact*sqrt(dfloat(jk+1)*dfloat(jb+1)*dfloat(lk2+1))        6d4s21
          isum=isb-lk2                                                   6d4s21
          if(mod(isum,2).ne.0)isum=isum+1                                6d4s21
          isum=isum/2                                                    6d4s21
          if(mod(isum,2).ne.0)fact=-fact                                 6d4s21
          ibra=joffb(jjb)+ioffb(jjb)+1                                  6d4s21
          do jrk=0,nrootk-1                                             6d4s21
           icol=iket+jrk
           do jrb=0,nrootb-1                                            6d4s21
            iadt=itmom(ipass)+jrb+nrootb*jrk                            6d4s21
            irow=ibra+jrb
            tm(irow,icol,isto)=bc(iadt)*fact                            6d4s21
           end do                                                       6d4s21
          end do                                                        6d4s21
          if(ipass.gt.1)then                                            6d4s21
           fact=f6j(jk,jb,iq2,lb2,lk2,isb,1)                              6d4s21
           fact=fact*sqrt(dfloat(jk+1)*dfloat(jb+1)*dfloat(lb2+1))        6d4s21
           isum=isb-lb2                                                   6d4s21
           if(mod(isum,2).ne.0)isum=isum+1                                6d4s21
           isum=isum/2                                                    6d4s21
           if(mod(isum,2).ne.0)fact=-fact                                 6d4s21
           do jrk=0,nrootk-1                                             6d4s21
            icol=iket+jrk
            do jrb=0,nrootb-1                                            6d4s21
             iadt=itmom(ipass)+jrb+nrootb*jrk                            6d4s21
             irow=ibra+jrb
             tm(icol,irow,isto)=bc(iadt)*fact                           6d4s21
            end do                                                       6d4s21
           end do                                                        6d4s21
          else                                                          6d4s21
           fact=f6j(jk,jb,iq2,lb2,lk2,isb,1)                              6d4s21
           fact=fact*sqrt(dfloat(jk+1)*dfloat(jb+1)*dfloat(lb2+1))        6d4s21
           isum=isb-lb2                                                   6d4s21
           if(mod(isum,2).ne.0)isum=isum+1                                6d4s21
           isum=isum/2                                                    6d4s21
           if(mod(isum,2).ne.0)fact=-fact                                 6d4s21
           do jrk=0,nrootk-1                                             6d4s21
            icol=iket+jrk
            do jrb=0,nrootb-1                                            6d4s21
             iadt=itmom(ipass)+jrb+nrootb*jrk                            6d4s21
             irow=ibra+jrb
             tmb(icol,irow)=bc(iadt)*fact                               6d4s21
            end do                                                       6d4s21
           end do                                                        6d4s21
          end if                                                        6d4s21
         end do
        end do
       end if                                                           6d4s21
      end do
      return
      end
