c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine atomll(iwavedat,nspc,eige,iqne,tmee,ndime,             6d4s21
     $     eigo,iqno,tmoo,ndimo,tmeo,tmoe,lprint,bc,ibc)                11d14s22
      implicit real*8 (a-h,o-z)                                         6d4s21
      real*4 data(2)
      integer*2 iqne(4,*),iqno(4,*)                                     6d4s21
      integer*2 ipack2(2)                                               6d4s21
      logical lprint                                                    3d2s22
      character*6 labb,labk
      character*4 chi,clo                                               6d4s21
      character*2 ttype                                                 6d4s21
      equivalence (ipack2,data(2)),(data,pack8)                         6d4s21
      dimension iwavedat(nspc,*),eige(*),tmee(ndime,ndime,2),           6d4s21
     $     eigo(*),tmoo(ndimo,ndimo,2),tmeo(ndime,*),tmoe(ndimo,*)      6d4s21
      include "common.store"                                            6d4s21
      if(lprint)then                                                    3d2s22
       write(6,*)('Hi, my name is atomll')
      end if                                                            3d2s22
      c=137.0359895d0                                                   6d4s21
      ci=1d0/c                                                          6d4s21
      ll=ibcoff                                                         6d4s21
      nll=0                                                             6d4s21
c     c=137 bohr/au time
c     c=3d10 cm/s
c     (3d10 cm/s)/(137 bohr/au time)
c     bohr=0.529177249d-8 cm
      aaa=3d10*ci/0.529177249d-8                                        6d4s21
      do ie=1,ndime                                                     6d4s21
       je2=iqne(3,ie)                                                   6d4s21
       if(mod(je2,2).eq.0)then                                          6d4s21
        jjh=je2/2                                                       6d4s21
        write(chi,2)jjh                                                 6d4s21
    2   format(i2,'  ')                                                 6d6s21
       else                                                             6d4s21
        write(chi,3)je2                                                 6d4s21
    3   format(i2,'/2')                                                 6d4s21
       end if                                                           6d4s21
       gei=1d0/dfloat(je2+1)                                            6d4s21
       iwfe=iqne(1,ie)                                                  6d4s21
       do l=1,6                                                         6d4s21
        if(iwavedat(13+l,iwfe).ne.0)then                                6d4s21
         labk(l:l)=char(iwavedat(13+l,iwfe))                             6d4s21
        else                                                            6d4s21
         labk(l:l)=' '                                                   6d4s21
        end if                                                          6d4s21
       end do                                                           6d4s21
       do io=1,ndimo                                                    6d4s21
        jo2=iqno(3,io)                                                  6d4s21
        if(mod(jo2,2).eq.0)then                                          6d4s21
         jjh=jo2/2                                                       6d4s21
         write(clo,2)jjh                                                 6d4s21
        else                                                             6d4s21
         write(clo,3)jo2                                                6d4s21
        end if                                                           6d4s21
        iwfo=iqno(1,io)                                                  6d4s21
        do l=1,6                                                         6d4s21
         if(iwavedat(13+l,iwfo).ne.0)then                                6d4s21
          labb(l:l)=char(iwavedat(13+l,iwfo))                             6d4s21
         else                                                            6d4s21
          labb(l:l)=' '                                                   6d4s21
         end if                                                          6d4s21
        end do                                                           6d4s21
        goi=1d0/dfloat(jo2+1)                                           6d4s21
        if(eige(ie).gt.eigo(io))then
         we=eige(ie)-eigo(io)
         wecm=we*219474.635d0
         if(we.ne.0d0)then                                              1d19s23
          wl=1d7/wecm
         else                                                           1d19s23
          wl=1d7                                                        1d19s23
         end if                                                         1d19s23
         strength=eina(tmoe(io,ie),gei,we,ci,1,aaa)
         if(strength.ne.0d0)then                                         6d4s21
          alt=strength/aaa
          data(1)=strength                                              6d4s21
          ipack2(1)=-io                                                 6d4s21
          ipack2(2)=ie                                                  6d4s21
          ibcoff=ll+nll                                                 6d4s21
          call enough('atomll.  1',bc,ibc)
          bc(ll+nll)=pack8                                              6d4s21
          nll=nll+1                                                     6d4s21
         end if
        else
         we=eigo(io)-eige(ie)
         wecm=we*219474.635d0
         wl=1d7/wecm
         strength=eina(tmeo(ie,io),goi,we,ci,1,aaa)
         if(strength.ne.0d0)then                                         6d4s21
          alt=strength/aaa
          data(1)=strength                                              6d4s21
          ipack2(2)=-io                                                 6d4s21
          ipack2(1)=ie                                                  6d4s21
          ibcoff=ll+nll                                                 6d4s21
          call enough('atomll.  2',bc,ibc)
          bc(ll+nll)=pack8                                              6d4s21
          nll=nll+1
         end if
        end if
    1    format(i6,3es15.7,2x,i1,1x,i1,a6,a4,' <-',2x,i1,1x,i1,a6,a4,    6d4s21
     $        1x,a2,es15.7)                                                    6d4s21
       end do                                                           6d4s21
      end do                                                            6d4s21
      do i=1,ndime                                                      6d4s21
       je2=iqne(3,i)                                                    6d4s21
       if(mod(je2,2).eq.0)then                                          6d4s21
        jjh=je2/2                                                       6d4s21
        write(chi,2)jjh                                                 6d4s21
       else                                                             6d4s21
        write(chi,3)je2                                                 6d4s21
       end if                                                           6d4s21
       gi=1d0/dfloat(je2+1)                                             6d4s21
       iwf=iqne(1,i)                                                    6d4s21
       do l=1,6                                                         6d4s21
        if(iwavedat(13+l,iwf).ne.0)then                                 6d4s21
         labk(l:l)=char(iwavedat(13+l,iwf))                             6d4s21
        else                                                            6d4s21
         labk(l:l)=' '                                                   6d4s21
        end if                                                          6d4s21
       end do                                                           6d4s21
       do j=1,i-1                                                       6d4s21
        jo2=iqne(3,j)                                                   6d4s21
        if(mod(jo2,2).eq.0)then                                          6d4s21
         jjh=jo2/2                                                       6d4s21
         write(clo,2)jjh                                                 6d4s21
        else                                                             6d4s21
         write(clo,3)jo2                                                 6d4s21
        end if                                                           6d4s21
        jwf=iqne(1,j)                                                   6d4s21
        do l=1,6                                                         6d4s21
         if(iwavedat(13+l,jwf).ne.0)then                                6d4s21
          labb(l:l)=char(iwavedat(13+l,jwf))                             6d4s21
         else                                                            6d4s21
          labb(l:l)=' '                                                   6d4s21
         end if                                                          6d4s21
        end do                                                           6d4s21
        gj=1d0/dfloat(jo2+1)                                            6d4s21
        if(eige(j).gt.eige(i))then                                      6d4s21
         we=eige(j)-eige(i)
         wecm=we*219474.635d0
         wl=1d7/wecm
         do ipass=1,2                                                   6d4s21
          strength=eina(tmee(i,j,ipass),gj,we,ci,1,aaa)
          if(strength.ne.0d0)then                                         6d4s21
           if(ipass.eq.1)then                                           6d6s21
            ttype='m1'                                                  6d6s21
            phs=1d0                                                     6d6s21
           else                                                         6d6s21
            ttype='e2'                                                  6d6s21
            phs=-1d0                                                    6d6s21
           end if                                                       6d6s21
           alt=strength/aaa
           data(1)=strength*phs                                         6d6s21
           ipack2(1)=i                                                  6d6s21
           ipack2(2)=j                                                  6d4s21
           ibcoff=ll+nll                                                 6d4s21
           call enough('atomll.  3',bc,ibc)
           bc(ll+nll)=pack8                                              6d4s21
           nll=nll+1                                                     6d4s21
          end if
         end do                                                         6d6s21
        else                                                            6d4s21
         we=eige(i)-eige(j)
         wecm=we*219474.635d0
         wl=1d7/wecm
         do ipass=1,2                                                   6d4s21
          strength=eina(tmee(j,i,ipass),gi,we,ci,1,aaa)
          if(strength.ne.0d0)then                                         6d4s21
           if(ipass.eq.1)then                                           6d6s21
            ttype='m1'                                                  6d6s21
            phs=1d0                                                     6d6s21
           else                                                         6d6s21
            ttype='e2'                                                  6d6s21
            phs=-1d0                                                    6d6s21
           end if                                                       6d6s21
           alt=strength/aaa
           data(1)=strength*phs                                         6d6s21
           ipack2(1)=j                                                  6d6s21
           ipack2(2)=i                                                  6d4s21
           ibcoff=ll+nll                                                 6d4s21
           call enough('atomll.  4',bc,ibc)
           bc(ll+nll)=pack8                                              6d4s21
           nll=nll+1                                                     6d4s21
          end if
         end do                                                         6d6s21
        end if                                                          6d4s21
       end do                                                           6d4s21
      end do                                                            6d4s21
      do i=1,ndimo                                                      6d4s21
       je2=iqno(3,i)                                                    6d4s21
       if(mod(je2,2).eq.0)then                                          6d4s21
        jjh=je2/2                                                       6d4s21
        write(chi,2)jjh                                                 6d4s21
       else                                                             6d4s21
        write(chi,3)je2                                                 6d4s21
       end if                                                           6d4s21
       gi=1d0/dfloat(je2+1)                                             6d4s21
       iwf=iqno(1,i)                                                    6d4s21
       do l=1,6                                                         6d4s21
        if(iwavedat(13+l,iwf).ne.0)then                                 6d4s21
         labk(l:l)=char(iwavedat(13+l,iwf))                             6d4s21
        else                                                            6d4s21
         labk(l:l)=' '                                                   6d4s21
        end if                                                          6d4s21
       end do                                                           6d4s21
       do j=1,i-1                                                       6d4s21
        jo2=iqno(3,j)                                                   6d4s21
        if(mod(jo2,2).eq.0)then                                          6d4s21
         jjh=jo2/2                                                       6d4s21
         write(clo,2)jjh                                                 6d4s21
        else                                                             6d4s21
         write(clo,3)jo2                                                 6d4s21
        end if                                                           6d4s21
        jwf=iqno(1,j)                                                   6d4s21
        do l=1,6                                                         6d4s21
         if(iwavedat(13+l,jwf).ne.0)then                                6d4s21
          labb(l:l)=char(iwavedat(13+l,jwf))                             6d4s21
         else                                                            6d4s21
          labb(l:l)=' '                                                   6d4s21
         end if                                                          6d4s21
        end do                                                           6d4s21
        gj=1d0/dfloat(jo2+1)                                            6d4s21
        if(eigo(j).gt.eigo(i))then                                      6d4s21
         we=eigo(j)-eigo(i)
         wecm=we*219474.635d0
         wl=1d7/wecm
         do ipass=1,2                                                   6d4s21
          strength=eina(tmoo(i,j,ipass),gj,we,ci,1,aaa)
          if(strength.ne.0d0)then                                         6d4s21
           if(ipass.eq.1)then                                           6d6s21
            ttype='m1'                                                  6d6s21
            phs=1d0                                                     6d6s21
           else                                                         6d6s21
            ttype='e2'                                                  6d6s21
            phs=-1d0                                                    6d6s21
           end if                                                       6d6s21
           alt=strength/aaa
           data(1)=strength*phs                                         6d6s21
           ipack2(1)=-i                                                  6d6s21
           ipack2(2)=-j                                                  6d4s21
           ibcoff=ll+nll                                                 6d4s21
           call enough('atomll.  5',bc,ibc)
           bc(ll+nll)=pack8                                              6d4s21
           nll=nll+1                                                     6d4s21
          end if
         end do                                                         6d6s21
        else                                                            6d4s21
         we=eigo(i)-eigo(j)
         wecm=we*219474.635d0
         wl=1d7/wecm
         do ipass=1,2                                                   6d4s21
          strength=eina(tmoo(j,i,ipass),gi,we,ci,1,aaa)
          if(strength.ne.0d0)then                                         6d4s21
           if(ipass.eq.1)then                                           6d6s21
            ttype='m1'                                                  6d6s21
            phs=1d0                                                     6d6s21
           else                                                         6d6s21
            ttype='e2'                                                  6d6s21
            phs=-1d0                                                    6d6s21
           end if                                                       6d6s21
           alt=strength/aaa
           data(1)=strength*phs                                         6d6s21
           ipack2(1)=-j                                                  6d6s21
           ipack2(2)=-i                                                  6d4s21
           ibcoff=ll+nll                                                 6d4s21
           call enough('atomll.  6',bc,ibc)
           bc(ll+nll)=pack8                                              6d4s21
           nll=nll+1                                                     6d4s21
          end if
         end do                                                         6d6s21
        end if                                                          6d4s21
       end do                                                           6d4s21
      end do                                                            6d4s21
      if(lprint)write(6,*)('total number of lines = '),nll              3d2s22
      isort=ll+nll                                                      6d4s21
      isort2=isort+nll                                                  6d4s21
      ibcoff=isort2+nll                                                 6d4s21
      call enough('atomll.  7',bc,ibc)
      do i=0,nll-1                                                      6d4s21
       pack8=bc(ll+i)                                                   6d4s21
       ihi=ipack2(2)                                                    6d4s21
       if(ihi.gt.0)then                                                 6d4s21
        ehi=eige(ihi)                                                   6d4s21
       else                                                             6d4s21
        ehi=eigo(-ihi)                                                  6d4s21
       end if                                                           6d4s21
       ilo=ipack2(1)                                                    6d4s21
       if(ilo.gt.0)then                                                 6d4s21
        elo=eige(ilo)
       else
        elo=eigo(-ilo)                                                  6d4s21
       end if                                                           6d4s21
       de=ehi-elo                                                       6d4s21
       bc(isort+i)=-de
      end do
      call dsortdws(bc(isort),ibc(isort2),nll)                          1d18s23
      if(lprint)write(6,10)                                             3d2s22
   10 format(14x,'transition e/wl',10x,'Einstein A',7x,'lower',13x,     6d6s21
     $     'upper',/,                                                   6d6s21
     $     16x,'1/cm',13x,'nm',12x,'1/s',2x,'rt',9x,'j',6x,'rt',9x,'j',  6d6s21
     $     2x,'type',2x,'elow(1/cm)')                                      6d6s21
      do i=0,nll-1                                                      6d4s21
       ip=i+1                                                           12d1s22
       ii=ibc(isort2+i)-1                                               6d4s21
       pack8=bc(ll+ii)                                                  6d4s21
       ilo=ipack2(1)                                                    6d4s21
       ihi=ipack2(2)                                                    6d4s21
       if(ihi.gt.0)then                                                 6d4s21
        ehi=eige(ihi)                                                   6d4s21
        iwh=iqne(1,ihi)                                                 6d4s21
        irh=iqne(2,ihi)                                                 6d4s21
        jh=iqne(3,ihi)                                                  6d4s21
       else                                                             6d4s21
        ehi=eigo(-ihi)                                                  6d4s21
        iwh=iqno(1,-ihi)                                                 6d4s21
        irh=iqno(2,-ihi)                                                 6d4s21
        jh=iqno(3,-ihi)                                                  6d4s21
       end if                                                           6d4s21
       do l=1,6                                                         6d4s21
        if(iwavedat(13+l,iwh).ne.0)then                                 6d4s21
         labb(l:l)=char(iwavedat(13+l,iwh))                             6d4s21
        else                                                            6d4s21
         labb(l:l)=' '                                                   6d4s21
        end if                                                          6d4s21
       end do                                                           6d4s21
       if(mod(jh,2).eq.0)then
        jhh=jh/2
        write(chi,2)jhh
       else
        write(chi,3)jh                                                  11d17s22
       end if
       ilo=ipack2(1)                                                    6d4s21
       if(ilo.gt.0)then                                                 6d4s21
        elo=eige(ilo)
        iwl=iqne(1,ilo)                                                 6d4s21
        irl=iqne(2,ilo)                                                 6d4s21
        jl=iqne(3,ilo)                                                  6d4s21
       else
        elo=eigo(-ilo)                                                  6d4s21
        iwl=iqno(1,-ilo)                                                 6d4s21
        irl=iqno(2,-ilo)                                                 6d4s21
        jl=iqno(3,-ilo)                                                  6d4s21
       end if                                                           6d4s21
       do l=1,6                                                         6d4s21
        if(iwavedat(13+l,iwl).ne.0)then                                 6d4s21
         labk(l:l)=char(iwavedat(13+l,iwl))                             6d4s21
        else                                                            6d4s21
         labk(l:l)=' '                                                   6d4s21
        end if                                                          6d4s21
       end do                                                           6d4s21
       if(mod(jl,2).eq.0)then
        jhh=jl/2
        write(clo,2)jhh
       else
        write(clo,3)jl                                                  11d21s22
       end if
       de=ehi-elo                                                       6d4s21
       if(ilo*ihi.gt.0)then                                             6d4s21
        if(data(1).gt.0d0)then                                          6d4s21
         ttype='m1'                                                     6d4s21
        else                                                            6d4s21
         ttype='e2'                                                     6d4s21
         data(1)=-data(1)                                               6d4s21
        end if                                                          6d4s21
       else                                                             6d4s21
        ttype='e1'                                                      6d4s21
       end if                                                           6d4s21
       decm=de*219474.635d0
       elo=elo*219474.635d0                                             6d6s21
       wl=1d7/decm
       if(lprint)write(6,1)ip,decm,wl,data(1),irl,iwavedat(1,iwl),labk, 12d1s22
     $      clo,irh,iwavedat(1,iwh),labb,chi,ttype,elo                  3d2s22
      end do                                                            6d4s21
      ibcoff=ll
      return                                                            6d4s21
      end                                                               6d4s21
