c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine paraeri4c(natom,ngaus,ibdat,nbasis,fock,den,           2d24s20
     $                iblstor,ibsstor,nbb,idwsdeb,ascale,nbaslarge,     2d24s20
     $     nbassmall,ilsz,isnorm,nbtot,irtyp,bc,ibc)                    11d17s22
      implicit real*8 (a-h,o-z)
      external second
      integer*8 nn8,ipack8                                              2d24s20
      integer*2 ipack2(4)                                               2d24s20
      equivalence (ipack2,ipack8)                                       2d24s20
      include "common.hf"
      include "common.store"
      include "common.spher"
      logical log1                                                      6d6s18
      integer*8 iblstor,ibsstor,la8,lb8,lc8,ld8                         2d24s20
      include "common.basis"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      dimension fock(*),den(*),iblstor(*),ibsstor(*),irtyp(5)           2d25s20
      scaleb=1d0
      ibcoffo=ibcoff                                                    2d19s10
      irowg=1
      icolg=2
      igoal1=irowg+nbtot*(icolg-1)
      igoal2=icolg+nbtot*(irowg-1)
      pi=acos(-1d0)                                                     2d21s20
      nb2=nbtot*nbtot                                                   2d24s20
      oneover4pi=0.25d0/pi                                              2d21s20
c     ascale=1/(4*c*c), thus c = 0.5/sqrt(ascale)
c
      clight=0.5d0/sqrt(ascale)                                         2d21s20
      if(idwsdeb.gt.100)then                                            5d25s18
       write(6,*)('in paraeric '),nbaslarge,nbassmall,nbb,ibcoff
       write(6,*)('ibdat '),ibdat
      do i=0,ngaus-1
       ip=i+1
       write(6,*)ip,ibc(ibdat+i),(bc(ibdat+i+ngaus*j),j=1,2),
     $      bc(isnorm+i),
     $      ibc(ibdat+i+3*ngaus),(bc(ibdat+i+ngaus*j),j=4,6),
     $      (ibc(ibdat+i+ngaus*j),j=7,8)
      end do
       call prntm2(bc(ibdat),ngaus,8,ngaus)
      end if                                                            5d25s18
      call second(time1)
      itry=13*25
c
c     build shell pair order
c
      nnt=ngaus**4
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nnt                                                  2d24s20
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nnt*2                                                   2d19s10
      call enough('paraeri4c.  1',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus7=ngaus*7                                                    5d3s10
      ii12=0
      do i1=1,ngaus                                                     2d19s10
       ipack2(1)=i1
       do i2=1,ngaus
        i12sz=((ibc(jbdat+i1)+ibc(jbdat+i2)+2)**2)                      2d24s20
        ipack2(2)=i2                                                    2d24s20
        jj12=0                                                          2d24s20
        do j1=1,ngaus
         do j2=1,ngaus
           j12sz=((ibc(jbdat+j1)+ibc(jbdat+2)+2)**2)                    2d24s20
           ibc(j12)=-i12sz*j12sz
           ipack2(3)=j1                                                 2d24s20
           ipack2(4)=j2                                                 2d24s20
           ibc(j12+nnt)=ipack8                                          2d24s20
           j12=j12+1
          jj12=jj12+1
         end do
        end do
   50   continue
        ii12=ii12+1
       end do                                                           2d19s10
      end do                                                            2d19s10
      is12=i12+nnt*2                                                    2d24s20
      nn8=nnt                                                           2d24s20
      call idsortdws(ibc(i12),ibc(is12),nn8)                            1d18s23
      do i=0,nn8-1
       ipack8=ibc(i12+i+nnt)
      end do
      ii1=i12+nnt                                                        2d19s10
      jpair=ipair                                                       2d19s10
      do i=1,nnt                                                         2d19s10
       j=ibc(is12+i-1)-1                                                2d19s10
    1  format(i5,i8,2i5)                                                2d19s10
       ibc(jpair)=ibc(ii1+j)                                            2d24s20
       ipack8=ibc(ii1+j)
       jpair=jpair+1                                                    2d24s20
      end do                                                            2d19s10
      do i=1,nbb                                                        12d28s19
       fock(i)=0d0                                                       2d24s10
      end do                                                            2d19s10
      call second(time2)
      telap=time2-time1
      loopit=0                                                          5d25s18
      do i=1+mynowprog,nnt,mynprocg                                     2d24s20
       loopit=loopit+1                                                  5d25s18
       jpair=ipair+i-1                                                  2d24s20
       ipack8=ibc(jpair)                                                2d24s20
       ja=jbdat+ipack2(1)
       isnorma=isnorm+ipack2(1)-1
       ja2=ja+ngaus
       ja3=ja2+ngaus
       ja4=ja3+ngaus
       ja5=ja4+ngaus
       ja6=ja5+ngaus
       ja7=ja6+ngaus
       ja8=ja7+ngaus
       jb=jbdat+ipack2(2)
       isnormb=isnorm+ipack2(2)-1
       jb2=jb+ngaus
       jb3=jb2+ngaus
       jb4=jb3+ngaus
       jb5=jb4+ngaus
       jb6=jb5+ngaus
       jb7=jb6+ngaus
       jb8=jb7+ngaus
       jc=jbdat+ipack2(3)
       isnormc=isnorm+ipack2(3)-1
       jc2=jc+ngaus
       jc3=jc2+ngaus
       jc4=jc3+ngaus
       jc5=jc4+ngaus
       jc6=jc5+ngaus
       jc7=jc6+ngaus
       jc8=jc7+ngaus
       jd=jbdat+ipack2(4)
       isnormd=isnorm+ipack2(4)-1
       jd2=jd+ngaus
       jd3=jd2+ngaus
       jd4=jd3+ngaus
       jd5=jd4+ngaus
       jd6=jd5+ngaus
       jd7=jd6+ngaus
       jd8=jd7+ngaus
       nal=ibc(ilsz+ibc(ja))                                            2d24s20
       nbl=ibc(ilsz+ibc(jb))                                            2d24s20
       ncl=ibc(ilsz+ibc(jc))                                            2d24s20
       ndl=ibc(ilsz+ibc(jd))                                            2d24s20
       nas=ibc(ilsz+ibc(ja)+1)                                          2d24s20
       nbs=ibc(ilsz+ibc(jb)+1)                                          2d24s20
       ncs=ibc(ilsz+ibc(jc)+1)                                          2d24s20
       nds=ibc(ilsz+ibc(jd)+1)                                          2d24s20
       la8=ibc(ja)+1                                                    2d24s20
       lb8=ibc(jb)+1                                                    2d24s20
       lc8=ibc(jc)+1                                                    2d24s20
       ld8=ibc(jd)+1                                                    2d24s20
c
c     llll
c
       if(irtyp(1).ne.0)then                                            2d25s20
        nla=2*ibc(ja)+1                                                  2d24s20
        nlb=2*ibc(jb)+1                                                  2d24s20
        nlc=2*ibc(jc)+1                                                  2d24s20
        nld=2*ibc(jd)+1                                                  2d24s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
c
c     ab electron 1 cd electron 2
c
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d24s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        ihit=0
        jb1=ib1                                                          2d24s20
        do id=0,ndl-1                                                    2d24s20
         idp=id+iblstor(ipack2(4))                                       2d24s20
         idpp=idp+nbaslarge                                              2d24s20
         do ic=0,ncl-1                                                   2d24s20
          icp=ic+iblstor(ipack2(3))                                      2d24s20
          icpp=icp+nbaslarge                                             2d24s20
          jfockar=icp+1+nbtot*idp                                         2d24s20
          jfockai=jfockar+nb2                                            2d25s20
          jfockbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbi=jfockbr+nb2                                            2d24s20
          do ib=0,nbl-1                                                  2d24s20
           ibp=ib+iblstor(ipack2(2))                                     2d24s20
           ibpp=ibp+nbaslarge                                            2d24s20
           do ia=0,nal-1                                                 2d24s20
            iap=ia+iblstor(ipack2(1))                                    2d24s20
            iapp=iap+nbaslarge                                           2d24s20
            jdenar=iap+1+nbtot*ibp                                       2d25s20
            jdenai=jdenar+nb2                                              2d24s20
            jdenbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbi=jdenbr+nb2                                            2d24s20
            dsr=den(jdenar)+den(jdenbr)                                  2d25s20
            dsi=den(jdenai)+den(jdenai)                                  2d25s20
            termr=dsr*bc(jb1)                                            2d25s20
            termi=dsi*bc(jb1)                                            2d25s20
            fock(jfockar)=fock(jfockar)+termr                            2d25s20
            fock(jfockbr)=fock(jfockbr)+termr                            2d25s20
            fock(jfockai)=fock(jfockai)+termi                            2d25s20
            fock(jfockbi)=fock(jfockbi)+termi                            2d25s20
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            fock(kfockaar)=fock(kfockaar)-den(kdenaar)*bc(jb1)             2d25s20
            fock(kfockaai)=fock(kfockaai)-den(kdenaai)*bc(jb1)             2d25s20
            fock(kfockabr)=fock(kfockabr)-den(kdenbar)*bc(jb1)             2d25s20
            fock(kfockabi)=fock(kfockabi)-den(kdenbai)*bc(jb1)             2d25s20
            fock(kfockbar)=fock(kfockbar)-den(kdenabr)*bc(jb1)             2d25s20
            fock(kfockbai)=fock(kfockbai)-den(kdenabi)*bc(jb1)             2d25s20
            fock(kfockbbr)=fock(kfockbbr)-den(kdenbbr)*bc(jb1)             2d25s20
            fock(kfockbbi)=fock(kfockbbi)-den(kdenbbi)*bc(jb1)             2d25s20
            jb1=jb1+1                                                    2d24s20
           end do
          end do
         end do
        end do
        if(ihit.gt.0)call prntm2(bc(ib1),nal*nbl,ncl*ndl,nal*nbl)       2d26s20
       end if                                                           2d25s20
       if(irtyp(2).ne.0)then                                            2d25s20
c
c     ssll and llss
c
        nla=2*la8+1                                                     2d25s20
        nlb=2*lb8+1                                                     2d25s20
        nlc=2*ibc(jc)+1                                                  2d24s20
        nld=2*ibc(jd)+1                                                  2d24s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        jb1=ib1                                                          2d24s20
        do id=0,ndl-1                                                    2d24s20
         idp=id+iblstor(ipack2(4))                                       2d24s20
         idpp=idp+nbaslarge                                              2d24s20
         do ic=0,ncl-1                                                   2d24s20
          icp=ic+iblstor(ipack2(3))                                      2d24s20
          icpp=icp+nbaslarge                                             2d24s20
          jfockar=icp+1+nbtot*idp                                         2d24s20
          jfockai=jfockar+nb2                                            2d25s20
          jfockbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbi=jfockbr+nb2                                            2d24s20
          do ib=0,nbs-1                                                  2d24s20
           ibp=ib+ibsstor(ipack2(2))+2*nbaslarge                        2d25s20
           ibpp=ibp+nbassmall                                           2d25s20
           do ia=0,nas-1                                                 2d24s20
            iap=ia+ibsstor(ipack2(1))+2*nbaslarge                       2d25s20
            iapp=iap+nbassmall                                          2d25s20
            jdenar=iap+1+nbtot*ibp                                       2d25s20
            jdenai=jdenar+nb2                                              2d24s20
            jdenbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbi=jdenbr+nb2                                            2d24s20
            dsr=den(jdenar)+den(jdenbr)                                  2d25s20
            dsi=den(jdenai)+den(jdenai)                                  2d25s20
            termr=dsr*bc(jb1)                                            2d25s20
            termi=dsi*bc(jb1)                                            2d25s20
            fock(jfockar)=fock(jfockar)+termr                            2d25s20
            fock(jfockbr)=fock(jfockbr)+termr                            2d25s20
            fock(jfockai)=fock(jfockai)+termi                            2d25s20
            fock(jfockbi)=fock(jfockbi)+termi                            2d25s20
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            fock(kfockaar)=fock(kfockaar)-den(kdenaar)*bc(jb1)             2d25s20
            fock(kfockaai)=fock(kfockaai)-den(kdenaai)*bc(jb1)             2d25s20
            fock(kfockabr)=fock(kfockabr)-den(kdenbar)*bc(jb1)             2d25s20
            fock(kfockabi)=fock(kfockabi)-den(kdenbai)*bc(jb1)             2d25s20
            fock(kfockbar)=fock(kfockbar)-den(kdenabr)*bc(jb1)             2d25s20
            fock(kfockbai)=fock(kfockbai)-den(kdenabi)*bc(jb1)             2d25s20
            fock(kfockbbr)=fock(kfockbbr)-den(kdenbbr)*bc(jb1)             2d25s20
            fock(kfockbbi)=fock(kfockbbi)-den(kdenbbi)*bc(jb1)             2d25s20
            jb1=jb1+1                                                    2d24s20
           end do
          end do
         end do
        end do
        nla=2*ibc(ja)+1                                                 2d25s20
        nlb=2*ibc(jb)+1                                                 2d25s20
        nlc=2*lc8+1                                                     2d25s20
        nld=2*ld8+1                                                     2d25s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        jb1=ib1                                                          2d24s20
        do id=0,nds-1                                                    2d24s20
         idp=id+ibsstor(ipack2(4))+2*nbaslarge                          2d25s20
         idpp=idp+nbassmall                                             2d25s20
         do ic=0,ncs-1                                                   2d24s20
          icp=ic+ibsstor(ipack2(3))+2*nbaslarge                         2d25s20
          icpp=icp+nbassmall                                            2d25s20
          jfockar=icp+1+nbtot*idp                                         2d24s20
          jfockai=jfockar+nb2                                            2d25s20
          jfockbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbi=jfockbr+nb2                                            2d24s20
          do ib=0,nbl-1                                                  2d24s20
           ibp=ib+iblstor(ipack2(2))                                    2d25s20
           ibpp=ibp+nbaslarge                                           2d25s20
           do ia=0,nal-1                                                 2d24s20
            iap=ia+iblstor(ipack2(1))                                   2d25s20
            iapp=iap+nbaslarge                                          2d25s20
            jdenar=iap+1+nbtot*ibp                                       2d25s20
            jdenai=jdenar+nb2                                              2d24s20
            jdenbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbi=jdenbr+nb2                                            2d24s20
            dsr=den(jdenar)+den(jdenbr)                                  2d25s20
            dsi=den(jdenai)+den(jdenai)                                  2d25s20
            termr=dsr*bc(jb1)                                            2d25s20
            termi=dsi*bc(jb1)                                            2d25s20
            fock(jfockar)=fock(jfockar)+termr                            2d25s20
            fock(jfockbr)=fock(jfockbr)+termr                            2d25s20
            fock(jfockai)=fock(jfockai)+termi                            2d25s20
            fock(jfockbi)=fock(jfockbi)+termi                            2d25s20
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            fock(kfockaar)=fock(kfockaar)-den(kdenaar)*bc(jb1)             2d25s20
            fock(kfockaai)=fock(kfockaai)-den(kdenaai)*bc(jb1)             2d25s20
            fock(kfockabr)=fock(kfockabr)-den(kdenbar)*bc(jb1)             2d25s20
            fock(kfockabi)=fock(kfockabi)-den(kdenbai)*bc(jb1)             2d25s20
            fock(kfockbar)=fock(kfockbar)-den(kdenabr)*bc(jb1)             2d25s20
            fock(kfockbai)=fock(kfockbai)-den(kdenabi)*bc(jb1)             2d25s20
            fock(kfockbbr)=fock(kfockbbr)-den(kdenbbr)*bc(jb1)             2d25s20
            fock(kfockbbi)=fock(kfockbbi)-den(kdenbbi)*bc(jb1)             2d25s20
            jb1=jb1+1                                                    2d24s20
           end do
          end do
         end do
        end do
       end if                                                           2d25s20
c
c     ssss                                                              2d24s20
c
       if(irtyp(3).ne.0)then                                            2d25s20
        nla=2*la8+1                                                      2d24s20
        nlb=2*lb8+1                                                      2d24s20
        nlc=2*lc8+1                                                      2d24s20
        nld=2*ld8+1                                                      2d24s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d24s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        jb1=ib1                                                          2d24s20
        do id=0,nds-1                                                    2d24s20
         idp=id+ibsstor(ipack2(4))+2*nbaslarge                          2d25s20
         idpp=idp+nbassmall                                             2d25s20
         do ic=0,ncs-1                                                   2d24s20
          icp=ic+ibsstor(ipack2(3))+2*nbaslarge                         2d25s20
          icpp=icp+nbassmall                                            2d25s20
          jfockar=icp+1+nbtot*idp                                         2d24s20
          jfockai=jfockar+nb2                                            2d25s20
          jfockbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbi=jfockbr+nb2                                            2d24s20
          do ib=0,nbs-1                                                  2d24s20
           ibp=ib+ibsstor(ipack2(2))+2*nbaslarge                        2d25s20
           ibpp=ibp+nbassmall                                           2d25s20
           do ia=0,nas-1                                                 2d24s20
            iap=ia+ibsstor(ipack2(1))+2*nbaslarge                       2d25s20
            iapp=iap+nbassmall                                          2d25s20
            jdenar=iap+1+nbtot*ibp                                       2d25s20
            jdenai=jdenar+nb2                                              2d24s20
            jdenbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbi=jdenbr+nb2                                            2d24s20
            dsr=den(jdenar)+den(jdenbr)                                  2d25s20
            dsi=den(jdenai)+den(jdenai)                                  2d25s20
            termr=dsr*bc(jb1)                                            2d25s20
            termi=dsi*bc(jb1)                                            2d25s20
            fock(jfockar)=fock(jfockar)+termr                            2d25s20
            fock(jfockbr)=fock(jfockbr)+termr                            2d25s20
            fock(jfockai)=fock(jfockai)+termi                            2d25s20
            fock(jfockbi)=fock(jfockbi)+termi                            2d25s20
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            fock(kfockaar)=fock(kfockaar)-den(kdenaar)*bc(jb1)             2d25s20
            fock(kfockaai)=fock(kfockaai)-den(kdenaai)*bc(jb1)             2d25s20
            fock(kfockabr)=fock(kfockabr)-den(kdenbar)*bc(jb1)             2d25s20
            fock(kfockabi)=fock(kfockabi)-den(kdenbai)*bc(jb1)             2d25s20
            fock(kfockbar)=fock(kfockbar)-den(kdenabr)*bc(jb1)             2d25s20
            fock(kfockbai)=fock(kfockbai)-den(kdenabi)*bc(jb1)             2d25s20
            fock(kfockbbr)=fock(kfockbbr)-den(kdenbbr)*bc(jb1)             2d25s20
            fock(kfockbbi)=fock(kfockbbi)-den(kdenbbi)*bc(jb1)             2d25s20
            jb1=jb1+1                                                    2d24s20
           end do
          end do
         end do
        end do
       end if
       if(irtyp(4).ne.0)then                                            2d25s20
c
c     magnetic Breit term: alpha1 dot alpha2/rij.
c     there is no half here, for the retardation term yields an identical
c     contribution. Gaunt approximation multiplies this again by 2.
c
        nla=2*ibc(ja)+1                                                     2d25s20
        nlb=2*lb8+1                                                     2d25s20
        nlc=2*ibc(jc)+1                                                  2d24s20
        nld=2*ld8+1                                                     2d25s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $       *dfloat(irtyp(4))                                          2d25s20
        xnorm=-xnorm*scaleb
        if(irtyp(5).ne.0)xnorm=xnorm*0.5d0                              2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        jb1=ib1                                                          2d24s20
        do id=0,nds-1                                                    2d24s20
         idp=id+ibsstor(ipack2(4))+nbaslarge*2                          2d25s20
         idpp=idp+nbassmall                                             2d25s20
         do ic=0,ncl-1                                                   2d24s20
          icp=ic+iblstor(ipack2(3))                                      2d24s20
          icpp=icp+nbaslarge                                             2d24s20
          jfockaar=icp+1+nbtot*idp                                         2d24s20
          jfockaai=jfockaar+nb2                                            2d25s20
          jfockabr=icp+1+nbtot*idpp                                         2d24s20
          jfockabi=jfockabr+nb2                                            2d25s20
          jfockbar=icpp+1+nbtot*idp                                     2d25s20
          jfockbai=jfockbar+nb2                                            2d24s20
          jfockbbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbbi=jfockbbr+nb2                                            2d24s20
          do ib=0,nbs-1                                                  2d24s20
           ibp=ib+ibsstor(ipack2(2))+2*nbaslarge                        2d25s20
           ibpp=ibp+nbassmall                                           2d25s20
           do ia=0,nal-1                                                 2d24s20
            iap=ia+iblstor(ipack2(1))                                   2d25s20
            iapp=iap+nbaslarge                                          2d25s20
            jdenaar=iap+1+nbtot*ibp                                     2d25s20
            jdenaai=jdenaar+nb2                                         2d25s20
            jdenabr=iap+1+nbtot*ibpp                                    2d25s20
            jdenabi=jdenabr+nb2                                         2d25s20
            jdenbar=iapp+1+nbtot*ibp                                    2d25s20
            jdenbai=jdenbar+nb2                                         2d25s20
            jdenbbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbbi=jdenbbr+nb2                                            2d24s20
            fock(jfockabr)=fock(jfockabr)+2d0*den(jdenbar)*bc(jb1)      2d25s20
            fock(jfockabi)=fock(jfockabi)+2d0*den(jdenbai)*bc(jb1)      2d25s20
            fock(jfockbar)=fock(jfockbar)+2d0*den(jdenabr)*bc(jb1)      2d25s20
            fock(jfockbai)=fock(jfockbai)+2d0*den(jdenabi)*bc(jb1)      2d25s20
            termr=(den(jdenaar)-den(jdenbbr))*bc(jb1)                   2d25s20
            termi=(den(jdenaai)-den(jdenbbi))*bc(jb1)                   2d25s20
            fock(jfockaar)=fock(jfockaar)+termr                         2d25s20
            fock(jfockaai)=fock(jfockaai)+termi                         2d25s20
            fock(jfockbbr)=fock(jfockbbr)+termr                         2d26s20
            fock(jfockbbi)=fock(jfockbbi)+termi                         2d26s20
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            fock(kfockaar)=fock(kfockaar)                               2d26s20
     $           -(den(kdenaar)+2d0*den(kdenbbr))*bc(jb1)               2d26s20
            fock(kfockaai)=fock(kfockaai)                               2d26s20
     $           -(den(kdenaai)+2d0*den(kdenbbi))*bc(jb1)               2d26s20
            fock(kfockbbr)=fock(kfockbbr)                               2d26s20
     $           -(2d0*den(kdenaar)-den(kdenbbr))*bc(jb1)               2d26s20
            fock(kfockbbi)=fock(kfockbbi)                               2d26s20
     $           -(2d0*den(kdenaai)-den(kdenbbi))*bc(jb1)               2d26s20
            fock(kfockabr)=fock(kfockabr)-den(kdenbar)*bc(jb1)               2d25s20
            fock(kfockabi)=fock(kfockabi)-den(kdenbai)*bc(jb1)               2d25s20
            fock(kfockbar)=fock(kfockbar)-den(kdenabr)*bc(jb1)               2d25s20
            fock(kfockbai)=fock(kfockbai)-den(kdenabi)*bc(jb1)               2d25s20
            jb1=jb1+1                                                    2d24s20
           end do
          end do
         end do
        end do
        nla=2*ibc(ja)+1                                                     2d25s20
        nlb=2*lb8+1                                                     2d25s20
        nlc=2*lc8+1                                                     2d25s20
        nld=2*ibc(jd)+1                                                 2d25s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $       *dfloat(irtyp(4))                                          2d25s20
        xnorm=-xnorm*scaleb
        if(irtyp(5).ne.0)xnorm=xnorm*0.5d0                              2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        jb1=ib1                                                          2d24s20
        do id=0,ndl-1                                                    2d24s20
         idp=id+iblstor(ipack2(4))                                      2d25s20
         idpp=idp+nbaslarge                                             2d25s20
         do ic=0,ncs-1                                                   2d24s20
          icp=ic+ibsstor(ipack2(3))+2*nbaslarge                         2d25s20
          icpp=icp+nbassmall                                            2d25s20
          jfockaar=icp+1+nbtot*idp                                         2d24s20
          jfockaai=jfockaar+nb2                                            2d25s20
          jfockabr=icp+1+nbtot*idpp                                         2d24s20
          jfockabi=jfockabr+nb2                                            2d25s20
          jfockbar=icpp+1+nbtot*idp                                     2d25s20
          jfockbai=jfockbar+nb2                                            2d24s20
          jfockbbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbbi=jfockbbr+nb2                                            2d24s20
          do ib=0,nbs-1                                                  2d24s20
           ibp=ib+ibsstor(ipack2(2))+2*nbaslarge                        2d25s20
           ibpp=ibp+nbassmall                                           2d25s20
           do ia=0,nal-1                                                 2d24s20
            iap=ia+iblstor(ipack2(1))                                   2d25s20
            iapp=iap+nbaslarge                                          2d25s20
            jdenaar=iap+1+nbtot*ibp                                     2d25s20
            jdenaai=jdenaar+nb2                                         2d25s20
            jdenabr=iap+1+nbtot*ibpp                                    2d25s20
            jdenabi=jdenabr+nb2                                         2d25s20
            jdenbar=iapp+1+nbtot*ibp                                    2d25s20
            jdenbai=jdenbar+nb2                                         2d25s20
            jdenbbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbbi=jdenbbr+nb2                                            2d24s20
            fock(jfockabr)=fock(jfockabr)+2d0*den(jdenbar)*bc(jb1)      2d25s20
            fock(jfockabi)=fock(jfockabi)+2d0*den(jdenbai)*bc(jb1)      2d25s20
            fock(jfockbar)=fock(jfockbar)+2d0*den(jdenabr)*bc(jb1)      2d25s20
            fock(jfockbai)=fock(jfockbai)+2d0*den(jdenabi)*bc(jb1)      2d25s20
            termr=(den(jdenaar)-den(jdenbbr))*bc(jb1)                   2d25s20
            termi=(den(jdenaai)-den(jdenbbi))*bc(jb1)                   2d25s20
            fock(jfockaar)=fock(jfockaar)+termr                         2d25s20
            fock(jfockaai)=fock(jfockaai)+termi                         2d25s20
            fock(jfockbbr)=fock(jfockbbr)+termr                         2d26s20
            fock(jfockbbi)=fock(jfockbbi)+termi                         2d26s20
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            fock(kfockaar)=fock(kfockaar)                               2d26s20
     $           -(den(kdenaar)+2d0*den(kdenbbr))*bc(jb1)               2d26s20
            fock(kfockaai)=fock(kfockaai)                               2d26s20
     $           -(den(kdenaai)+2d0*den(kdenbbi))*bc(jb1)               2d26s20
            fock(kfockbbr)=fock(kfockbbr)                               2d26s20
     $           -(2d0*den(kdenaar)-den(kdenbbr))*bc(jb1)               2d26s20
            fock(kfockbbi)=fock(kfockbbi)                               2d26s20
     $           -(2d0*den(kdenaai)-den(kdenbbi))*bc(jb1)               2d26s20
            fock(kfockabr)=fock(kfockabr)-den(kdenbar)*bc(jb1)               2d25s20
            fock(kfockabi)=fock(kfockabi)-den(kdenbai)*bc(jb1)               2d25s20
            fock(kfockbar)=fock(kfockbar)-den(kdenabr)*bc(jb1)               2d25s20
            fock(kfockbai)=fock(kfockbai)-den(kdenabi)*bc(jb1)               2d25s20
            jb1=jb1+1                                                    2d24s20
           end do
          end do
         end do
        end do
        nla=2*la8+1                                                     2d25s20
        nlb=2*ibc(jb)+1                                                     2d25s20
        nlc=2*ibc(jc)+1                                                  2d24s20
        nld=2*ld8+1                                                     2d25s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $       *dfloat(irtyp(4))                                          2d25s20
        xnorm=-xnorm*scaleb
        if(irtyp(5).ne.0)xnorm=xnorm*0.5d0                              2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        jb1=ib1                                                          2d24s20
        do id=0,nds-1                                                    2d24s20
         idp=id+ibsstor(ipack2(4))+nbaslarge*2                          2d25s20
         idpp=idp+nbassmall                                             2d25s20
         do ic=0,ncl-1                                                   2d24s20
          icp=ic+iblstor(ipack2(3))                                      2d24s20
          icpp=icp+nbaslarge                                             2d24s20
          jfockaar=icp+1+nbtot*idp                                         2d24s20
          jfockaai=jfockaar+nb2                                            2d25s20
          jfockabr=icp+1+nbtot*idpp                                         2d24s20
          jfockabi=jfockabr+nb2                                            2d25s20
          jfockbar=icpp+1+nbtot*idp                                     2d25s20
          jfockbai=jfockbar+nb2                                            2d24s20
          jfockbbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbbi=jfockbbr+nb2                                            2d24s20
          do ib=0,nbl-1                                                  2d24s20
           ibp=ib+iblstor(ipack2(2))                                    2d25s20
           ibpp=ibp+nbaslarge                                           2d25s20
           do ia=0,nas-1                                                 2d24s20
            iap=ia+ibsstor(ipack2(1))+nbaslarge*2                       2d25s20
            iapp=iap+nbassmall                                          2d25s20
            jdenaar=iap+1+nbtot*ibp                                     2d25s20
            jdenaai=jdenaar+nb2                                         2d25s20
            jdenabr=iap+1+nbtot*ibpp                                    2d25s20
            jdenabi=jdenabr+nb2                                         2d25s20
            jdenbar=iapp+1+nbtot*ibp                                    2d25s20
            jdenbai=jdenbar+nb2                                         2d25s20
            jdenbbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbbi=jdenbbr+nb2                                            2d24s20
            fock(jfockabr)=fock(jfockabr)+2d0*den(jdenbar)*bc(jb1)      2d25s20
            fock(jfockabi)=fock(jfockabi)+2d0*den(jdenbai)*bc(jb1)      2d25s20
            fock(jfockbar)=fock(jfockbar)+2d0*den(jdenabr)*bc(jb1)      2d25s20
            fock(jfockbai)=fock(jfockbai)+2d0*den(jdenabi)*bc(jb1)      2d25s20
            termr=(den(jdenaar)-den(jdenbbr))*bc(jb1)                   2d25s20
            termi=(den(jdenaai)-den(jdenbbi))*bc(jb1)                   2d25s20
            fock(jfockaar)=fock(jfockaar)+termr                         2d25s20
            fock(jfockaai)=fock(jfockaai)+termi                         2d25s20
            fock(jfockbbr)=fock(jfockbbr)+termr                         2d26s20
            fock(jfockbbi)=fock(jfockbbi)+termi                         2d26s20
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            fock(kfockaar)=fock(kfockaar)                               2d26s20
     $           -(den(kdenaar)+2d0*den(kdenbbr))*bc(jb1)               2d26s20
            fock(kfockaai)=fock(kfockaai)                               2d26s20
     $           -(den(kdenaai)+2d0*den(kdenbbi))*bc(jb1)               2d26s20
            fock(kfockbbr)=fock(kfockbbr)                               2d26s20
     $           -(2d0*den(kdenaar)-den(kdenbbr))*bc(jb1)               2d26s20
            fock(kfockbbi)=fock(kfockbbi)                               2d26s20
     $           -(2d0*den(kdenaai)-den(kdenbbi))*bc(jb1)               2d26s20
            fock(kfockabr)=fock(kfockabr)-den(kdenbar)*bc(jb1)               2d25s20
            fock(kfockabi)=fock(kfockabi)-den(kdenbai)*bc(jb1)               2d25s20
            fock(kfockbar)=fock(kfockbar)-den(kdenabr)*bc(jb1)               2d25s20
            fock(kfockbai)=fock(kfockbai)-den(kdenabi)*bc(jb1)               2d25s20
            jb1=jb1+1                                                    2d24s20
           end do
          end do
         end do
        end do
        nla=2*la8+1                                                     2d25s20
        nlb=2*ibc(jb)+1                                                     2d25s20
        nlc=2*lc8+1                                                     2d25s20
        nld=2*ibc(jd)+1                                                 2d25s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $       *dfloat(irtyp(4))                                          2d25s20
        xnorm=-xnorm*scaleb
        if(irtyp(5).ne.0)xnorm=xnorm*0.5d0                              2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d25s20
     $      dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        jb1=ib1                                                          2d24s20
        do id=0,ndl-1                                                    2d24s20
         idp=id+iblstor(ipack2(4))                                      2d25s20
         idpp=idp+nbaslarge                                             2d25s20
         do ic=0,ncs-1                                                   2d24s20
          icp=ic+ibsstor(ipack2(3))+2*nbaslarge                         2d25s20
          icpp=icp+nbassmall                                            2d25s20
          jfockaar=icp+1+nbtot*idp                                         2d24s20
          jfockaai=jfockaar+nb2                                            2d25s20
          jfockabr=icp+1+nbtot*idpp                                         2d24s20
          jfockabi=jfockabr+nb2                                            2d25s20
          jfockbar=icpp+1+nbtot*idp                                     2d25s20
          jfockbai=jfockbar+nb2                                            2d24s20
          jfockbbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbbi=jfockbbr+nb2                                            2d24s20
          do ib=0,nbl-1                                                  2d24s20
           ibp=ib+iblstor(ipack2(2))                                    2d25s20
           ibpp=ibp+nbaslarge                                           2d25s20
           do ia=0,nas-1                                                 2d24s20
            iap=ia+ibsstor(ipack2(1))+nbaslarge*2                       2d25s20
            iapp=iap+nbassmall                                          2d25s20
            jdenaar=iap+1+nbtot*ibp                                     2d25s20
            jdenaai=jdenaar+nb2                                         2d25s20
            jdenabr=iap+1+nbtot*ibpp                                    2d25s20
            jdenabi=jdenabr+nb2                                         2d25s20
            jdenbar=iapp+1+nbtot*ibp                                    2d25s20
            jdenbai=jdenbar+nb2                                         2d25s20
            jdenbbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbbi=jdenbbr+nb2                                            2d24s20
            fock(jfockabr)=fock(jfockabr)+2d0*den(jdenbar)*bc(jb1)      2d25s20
            fock(jfockabi)=fock(jfockabi)+2d0*den(jdenbai)*bc(jb1)      2d25s20
            fock(jfockbar)=fock(jfockbar)+2d0*den(jdenabr)*bc(jb1)      2d25s20
            fock(jfockbai)=fock(jfockbai)+2d0*den(jdenabi)*bc(jb1)      2d25s20
            termr=(den(jdenaar)-den(jdenbbr))*bc(jb1)                   2d25s20
            termi=(den(jdenaai)-den(jdenbbi))*bc(jb1)                   2d25s20
            fock(jfockaar)=fock(jfockaar)+termr                         2d25s20
            fock(jfockaai)=fock(jfockaai)+termi                         2d25s20
            fock(jfockbbr)=fock(jfockbbr)+termr                         2d26s20
            fock(jfockbbi)=fock(jfockbbi)+termi                         2d26s20
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            fock(kfockaar)=fock(kfockaar)                               2d26s20
     $           -(den(kdenaar)+2d0*den(kdenbbr))*bc(jb1)               2d26s20
            fock(kfockaai)=fock(kfockaai)                               2d26s20
     $           -(den(kdenaai)+2d0*den(kdenbbi))*bc(jb1)               2d26s20
            fock(kfockbbr)=fock(kfockbbr)                               2d26s20
     $           -(2d0*den(kdenaar)-den(kdenbbr))*bc(jb1)               2d26s20
            fock(kfockbbi)=fock(kfockbbi)                               2d26s20
     $           -(2d0*den(kdenaai)-den(kdenbbi))*bc(jb1)               2d26s20
            fock(kfockabr)=fock(kfockabr)-den(kdenbar)*bc(jb1)               2d25s20
            fock(kfockabi)=fock(kfockabi)-den(kdenbai)*bc(jb1)               2d25s20
            fock(kfockbar)=fock(kfockbar)-den(kdenabr)*bc(jb1)               2d25s20
            fock(kfockbai)=fock(kfockbai)-den(kdenabi)*bc(jb1)               2d25s20
            jb1=jb1+1                                                    2d24s20
           end do
          end do
         end do
        end do
       end if                                                           2d25s20
       if(irtyp(5).ne.0)then                                            2d26s20
c
c     retardation term.
c
        nla=2*ibc(ja)+1                                                     2d25s20
        nlb=2*lb8+1                                                     2d25s20
        nlc=2*ibc(jc)+1                                                  2d24s20
        nld=2*ld8+1                                                     2d25s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $       *dfloat(irtyp(5))                                          2d26s20
        xnorm=-xnorm*0.5d0*scaleb                                              2d26s20
        nabcd=nal*nbs*ncl*nds                                           2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,inorm,0,0,0,0,0,0,0,0,0,0,0,0,                   2d26s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=inorm+nabcd                                              2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,izz1,0,0,1,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=izz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,izz2,0,0,0,0,0,1,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(izz1+ii)=bc(izz1+ii)+bc(izz2+ii)+bc(inorm+ii)               2d26s20
        end do                                                          2d26s20
        ibcoff=izz2                                                     2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixx1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         2,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixx1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixx2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         2,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixx2+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyy1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=iyy1+nabcd
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyy2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixx1+ii)=bc(ixx1+ii)+bc(ixx2+ii)+bc(iyy1+ii)+bc(iyy2+ii)    2d26s20
     $        +2d0*bc(inorm+ii)                                         2d26s20
        end do
        ibcoff=ixx2
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixz1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixz2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixz1+ii)=bc(ixz1+ii)+bc(ixz2+ii)                               2d26s20
        end do                                                          2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyz1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=iyz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyz2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(iyz1+ii)=bc(iyz1+ii)+bc(iyz2+ii)                               2d26s20
        end do                                                          2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixy1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixy1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixy2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixy1+ii)=bc(ixy1+ii)+bc(ixy2+ii)                               2d26s20
        end do                                                          2d26s20
        jzz1=izz1                                                       2d26s20
        jxx1=ixx1                                                       2d26s20
        jxz1=ixz1                                                       2d26s20
        jyz1=iyz1                                                       2d26s20
        jxy1=ixy1                                                       2d26s20
        do id=0,nds-1                                                    2d24s20
         idp=id+ibsstor(ipack2(4))+nbaslarge*2                          2d25s20
         idpp=idp+nbassmall                                             2d25s20
         do ic=0,ncl-1                                                   2d24s20
          icp=ic+iblstor(ipack2(3))                                      2d24s20
          icpp=icp+nbaslarge                                             2d24s20
          jfockaar=icp+1+nbtot*idp                                         2d24s20
          jfockaai=jfockaar+nb2                                            2d25s20
          jfockabr=icp+1+nbtot*idpp                                         2d24s20
          jfockabi=jfockabr+nb2                                            2d25s20
          jfockbar=icpp+1+nbtot*idp                                     2d25s20
          jfockbai=jfockbar+nb2                                            2d24s20
          jfockbbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbbi=jfockbbr+nb2                                            2d24s20
          do ib=0,nbs-1                                                  2d24s20
           ibp=ib+ibsstor(ipack2(2))+2*nbaslarge                        2d25s20
           ibpp=ibp+nbassmall                                           2d25s20
           do ia=0,nal-1                                                 2d24s20
            iap=ia+iblstor(ipack2(1))                                   2d25s20
            iapp=iap+nbaslarge                                          2d25s20
            jdenaar=iap+1+nbtot*ibp                                     2d25s20
            jdenaai=jdenaar+nb2                                         2d25s20
            jdenabr=iap+1+nbtot*ibpp                                    2d25s20
            jdenabi=jdenabr+nb2                                         2d25s20
            jdenbar=iapp+1+nbtot*ibp                                    2d25s20
            jdenbai=jdenbar+nb2                                         2d25s20
            jdenbbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbbi=jdenbbr+nb2                                            2d24s20
            termr=bc(jzz1)*(den(jdenaar)-den(jdenbbr))                  2d26s20
     $           +bc(jxz1)*(den(jdenabr)+den(jdenbar))                  2d26s20
     $           +bc(jyz1)*(-den(jdenbai)+den(jdenabi))                 2d26s20
            termi=bc(jzz1)*(den(jdenaai)-den(jdenbbi))                  2d26s20
     $           +bc(jxz1)*(den(jdenabi)+den(jdenbai))                  2d26s20
     $           +bc(jyz1)*(den(jdenbar)-den(jdenabr))                  2d26s20
            fock(jfockaar)=fock(jfockaar)+termr                         2d25s20
            fock(jfockaai)=fock(jfockaai)+termi                         2d25s20
            fock(jfockbbr)=fock(jfockbbr)-termr                         2d26s20
            fock(jfockbbi)=fock(jfockbbi)-termi                         2d26s20
            ddr=den(jdenaar)-den(jdenbbr)                               2d26s20
            ddi=den(jdenaai)-den(jdenbbi)                               2d26s20
            termr1=bc(jxz1)*ddr+bc(jxx1)*(den(jdenbar)+den(jdenabr))
            termr2=-bc(jyz1)*ddi                                        2d26s20
            termi1=bc(jxz1)*ddi+bc(jxx1)*(den(jdenbai)+den(jdenabi))
            termi2=bc(jyz1)*ddr
            fock(jfockabr)=fock(jfockabr)+termr1-termr2
     $           +2d0*bc(jxy1)*den(jdenabi)
            fock(jfockabi)=fock(jfockabi)+termi1-termi2
     $           -2d0*bc(jxy1)*den(jdenabr)
            fock(jfockbar)=fock(jfockbar)+termr1+termr2
     $           -2d0*bc(jxy1)*den(jdenbai)
            fock(jfockbai)=fock(jfockbai)+termi1+termi2
     $           +2d0*bc(jxy1)*den(jdenbar)
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            termr=bc(jxz1)*(den(kdenbar)+den(kdenabr))
     $           -bc(jyz1)*(den(kdenbai)-den(kdenabi))
            termi=bc(jxz1)*(den(kdenbai)+den(kdenabi))                  2d26s20
     $           +bc(jyz1)*(den(kdenbar)-den(kdenabr))                  2d26s20
            fock(kfockaar)=fock(kfockaar)-termr-den(kdenaar)*bc(jzz1)   2d26s20
     $           -den(kdenbbr)*bc(jxx1)                                 2d26s20
            fock(kfockaai)=fock(kfockaai)-termi-den(kdenaai)*bc(jzz1)   2d26s20
     $           -den(kdenbbi)*bc(jxx1)                                 2d26s20
            fock(kfockbbr)=fock(kfockbbr)+termr-den(kdenaar)*bc(jxx1)   2d26s20
     $           -den(kdenbbr)*bc(jzz1)                                 2d26s20
            fock(kfockbbi)=fock(kfockbbi)+termi-den(kdenaai)*bc(jxx1)   2d26s20
     $           -den(kdenbbi)*bc(jzz1)                                 2d26s20
            ddr=den(kdenaar)-den(kdenbbr)                               2d26s20
            ddi=den(kdenaai)-den(kdenbbi)                               2d26s20
            termr1=bc(jxz1)*ddr                                         2d26s20
            termr2=-bc(jyz1)*ddi                                        2d26s20
            termi1=bc(jxz1)*ddi                                         2d26s20
            termi2=bc(jyz1)*ddr
            fock(kfockabr)=fock(kfockabr)-termr1+termr2
     $           -bc(jxx1)*den(kdenabr)+bc(jzz1)*den(kdenbar)
     $           -2d0*bc(jxy1)*den(kdenabi)                             2d26s20
            fock(kfockabi)=fock(kfockabi)-termi1+termi2
     $           -bc(jxx1)*den(kdenabi)+bc(jzz1)*den(kdenbai)
     $           +2d0*bc(jxy1)*den(kdenabr)                             2d26s20
            fock(kfockbar)=fock(kfockbar)-termr1-termr2
     $           -bc(jxx1)*den(kdenbar)+bc(jzz1)*den(kdenabr)
     $           +2d0*bc(jxy1)*den(kdenbai)                             2d26s20
            fock(kfockbai)=fock(kfockbai)-termi1-termi2
     $           -bc(jxx1)*den(kdenbai)+bc(jzz1)*den(kdenabi)
     $           -2d0*bc(jxy1)*den(kdenbar)                             2d26s20
            jzz1=jzz1+1                                                    2d24s20
            jxx1=jxx1+1                                                    2d24s20
            jxz1=jxz1+1                                                    2d24s20
            jyz1=jyz1+1                                                    2d24s20
            jxy1=jxy1+1                                                    2d24s20
           end do
          end do
         end do
        end do
        ibcoff=inorm                                                    2d26s20
        nla=2*ibc(ja)+1                                                     2d25s20
        nlb=2*lb8+1                                                     2d25s20
        nlc=2*lc8+1                                                     2d26s20
        nld=2*ibc(jd)+1                                                 2d26s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $       *dfloat(irtyp(5))                                          2d26s20
        xnorm=-xnorm*0.5d0*scaleb                                              2d26s20
        nabcd=nal*nbs*ncl*nds                                           2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,inorm,0,0,0,0,0,0,0,0,0,0,0,0,                    2d26s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=inorm+nabcd                                              2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,izz1,0,0,1,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=izz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,izz2,0,0,0,0,0,1,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(izz1+ii)=bc(izz1+ii)+bc(izz2+ii)+bc(inorm+ii)               2d26s20
        end do                                                          2d26s20
        ibcoff=izz2                                                     2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixx1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         2,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixx1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixx2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         2,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixx2+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyy1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=iyy1+nabcd
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyy2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixx1+ii)=bc(ixx1+ii)+bc(ixx2+ii)+bc(iyy1+ii)+bc(iyy2+ii)    2d26s20
     $        +2d0*bc(inorm+ii)                                         2d26s20
        end do
        ibcoff=ixx2
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixz1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixz2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixz1+ii)=bc(ixz1+ii)+bc(ixz2+ii)                               2d26s20
        end do                                                          2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyz1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=iyz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyz2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(iyz1+ii)=bc(iyz1+ii)+bc(iyz2+ii)                               2d26s20
        end do                                                          2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixy1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixy1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixy2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixy1+ii)=bc(ixy1+ii)+bc(ixy2+ii)                               2d26s20
        end do                                                          2d26s20
        jzz1=izz1                                                       2d26s20
        jxx1=ixx1                                                       2d26s20
        jxz1=ixz1                                                       2d26s20
        jyz1=iyz1                                                       2d26s20
        jxy1=ixy1                                                       2d26s20
        do id=0,ndl-1                                                    2d24s20
         idp=id+iblstor(ipack2(4))                                      2d26s20
         idpp=idp+nbaslarge
         do ic=0,ncs-1                                                   2d24s20
          icp=ic+ibsstor(ipack2(3))+2*nbaslarge                         2d26s20
          icpp=icp+nbassmall                                            2d26s20
          jfockaar=icp+1+nbtot*idp                                         2d24s20
          jfockaai=jfockaar+nb2                                            2d25s20
          jfockabr=icp+1+nbtot*idpp                                         2d24s20
          jfockabi=jfockabr+nb2                                            2d25s20
          jfockbar=icpp+1+nbtot*idp                                     2d25s20
          jfockbai=jfockbar+nb2                                            2d24s20
          jfockbbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbbi=jfockbbr+nb2                                            2d24s20
          do ib=0,nbs-1                                                  2d24s20
           ibp=ib+ibsstor(ipack2(2))+2*nbaslarge                        2d25s20
           ibpp=ibp+nbassmall                                           2d25s20
           do ia=0,nal-1                                                 2d24s20
            iap=ia+iblstor(ipack2(1))                                   2d25s20
            iapp=iap+nbaslarge                                          2d25s20
            jdenaar=iap+1+nbtot*ibp                                     2d25s20
            jdenaai=jdenaar+nb2                                         2d25s20
            jdenabr=iap+1+nbtot*ibpp                                    2d25s20
            jdenabi=jdenabr+nb2                                         2d25s20
            jdenbar=iapp+1+nbtot*ibp                                    2d25s20
            jdenbai=jdenbar+nb2                                         2d25s20
            jdenbbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbbi=jdenbbr+nb2                                            2d24s20
            termr=bc(jzz1)*(den(jdenaar)-den(jdenbbr))                  2d26s20
     $           +bc(jxz1)*(den(jdenabr)+den(jdenbar))                  2d26s20
     $           +bc(jyz1)*(-den(jdenbai)+den(jdenabi))                 2d26s20
            termi=bc(jzz1)*(den(jdenaai)-den(jdenbbi))                  2d26s20
     $           +bc(jxz1)*(den(jdenabi)+den(jdenbai))                  2d26s20
     $           +bc(jyz1)*(den(jdenbar)-den(jdenabr))                  2d26s20
            fock(jfockaar)=fock(jfockaar)+termr                         2d25s20
            fock(jfockaai)=fock(jfockaai)+termi                         2d25s20
            fock(jfockbbr)=fock(jfockbbr)-termr                         2d26s20
            fock(jfockbbi)=fock(jfockbbi)-termi                         2d26s20
            ddr=den(jdenaar)-den(jdenbbr)                               2d26s20
            ddi=den(jdenaai)-den(jdenbbi)                               2d26s20
            termr1=bc(jxz1)*ddr+bc(jxx1)*(den(jdenbar)+den(jdenabr))
            termr2=-bc(jyz1)*ddi                                        2d26s20
            termi1=bc(jxz1)*ddi+bc(jxx1)*(den(jdenbai)+den(jdenabi))
            termi2=bc(jyz1)*ddr
            fock(jfockabr)=fock(jfockabr)+termr1-termr2
     $           +2d0*bc(jxy1)*den(jdenabi)
            fock(jfockabi)=fock(jfockabi)+termi1-termi2
     $           -2d0*bc(jxy1)*den(jdenabr)
            fock(jfockbar)=fock(jfockbar)+termr1+termr2
     $           -2d0*bc(jxy1)*den(jdenbai)
            fock(jfockbai)=fock(jfockbai)+termi1+termi2
     $           +2d0*bc(jxy1)*den(jdenbar)
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            termr=bc(jxz1)*(den(kdenbar)+den(kdenabr))
     $           -bc(jyz1)*(den(kdenbai)-den(kdenabi))
            termi=bc(jxz1)*(den(kdenbai)+den(kdenabi))                  2d26s20
     $           +bc(jyz1)*(den(kdenbar)-den(kdenabr))                  2d26s20
            fock(kfockaar)=fock(kfockaar)-termr-den(kdenaar)*bc(jzz1)   2d26s20
     $           -den(kdenbbr)*bc(jxx1)                                 2d26s20
            fock(kfockaai)=fock(kfockaai)-termi-den(kdenaai)*bc(jzz1)   2d26s20
     $           -den(kdenbbi)*bc(jxx1)                                 2d26s20
            fock(kfockbbr)=fock(kfockbbr)+termr-den(kdenaar)*bc(jxx1)   2d26s20
     $           -den(kdenbbr)*bc(jzz1)                                 2d26s20
            fock(kfockbbi)=fock(kfockbbi)+termi-den(kdenaai)*bc(jxx1)   2d26s20
     $           -den(kdenbbi)*bc(jzz1)                                 2d26s20
            ddr=den(kdenaar)-den(kdenbbr)                               2d26s20
            ddi=den(kdenaai)-den(kdenbbi)                               2d26s20
            termr1=bc(jxz1)*ddr                                         2d26s20
            termr2=-bc(jyz1)*ddi                                        2d26s20
            termi1=bc(jxz1)*ddi                                         2d26s20
            termi2=bc(jyz1)*ddr
            fock(kfockabr)=fock(kfockabr)-termr1+termr2
     $           -bc(jxx1)*den(kdenabr)+bc(jzz1)*den(kdenbar)
     $           -2d0*bc(jxy1)*den(kdenabi)                             2d26s20
            fock(kfockabi)=fock(kfockabi)-termi1+termi2
     $           -bc(jxx1)*den(kdenabi)+bc(jzz1)*den(kdenbai)
     $           +2d0*bc(jxy1)*den(kdenabr)                             2d26s20
            fock(kfockbar)=fock(kfockbar)-termr1-termr2
     $           -bc(jxx1)*den(kdenbar)+bc(jzz1)*den(kdenabr)
     $           +2d0*bc(jxy1)*den(kdenbai)                             2d26s20
            fock(kfockbai)=fock(kfockbai)-termi1-termi2
     $           -bc(jxx1)*den(kdenbai)+bc(jzz1)*den(kdenabi)
     $           -2d0*bc(jxy1)*den(kdenbar)                             2d26s20
            jzz1=jzz1+1                                                    2d24s20
            jxx1=jxx1+1                                                    2d24s20
            jxz1=jxz1+1                                                    2d24s20
            jyz1=jyz1+1                                                    2d24s20
            jxy1=jxy1+1                                                    2d24s20
           end do
          end do
         end do
        end do
        ibcoff=inorm                                                    2d26s20
        nla=2*la8+1                                                     2d25s20
        nlb=2*ibc(jb)+1                                                 2d26s20
        nlc=2*ibc(jc)+1                                                  2d24s20
        nld=2*ld8+1                                                     2d25s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $       *dfloat(irtyp(5))                                          2d26s20
        xnorm=-xnorm*0.5d0*scaleb                                              2d26s20
        nabcd=nal*nbs*ncl*nds                                           2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,inorm,0,0,0,0,0,0,0,0,0,0,0,0,                   2d26s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=inorm+nabcd                                              2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,izz1,0,0,1,0,0,0,0,0,0,0,0,0,                    2d26s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=izz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,izz2,0,0,0,0,0,1,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(izz1+ii)=bc(izz1+ii)+bc(izz2+ii)+bc(inorm+ii)               2d26s20
        end do                                                          2d26s20
        ibcoff=izz2                                                     2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixx1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         2,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixx1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixx2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         2,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixx2+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyy1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=iyy1+nabcd
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyy2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixx1+ii)=bc(ixx1+ii)+bc(ixx2+ii)+bc(iyy1+ii)+bc(iyy2+ii)    2d26s20
     $        +2d0*bc(inorm+ii)                                         2d26s20
        end do
        ibcoff=ixx2
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixz1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixz2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixz1+ii)=bc(ixz1+ii)+bc(ixz2+ii)                               2d26s20
        end do                                                          2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyz1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=iyz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyz2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(iyz1+ii)=bc(iyz1+ii)+bc(iyz2+ii)                               2d26s20
        end do                                                          2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixy1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixy1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixy2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixy1+ii)=bc(ixy1+ii)+bc(ixy2+ii)                               2d26s20
        end do                                                          2d26s20
        jzz1=izz1                                                       2d26s20
        jxx1=ixx1                                                       2d26s20
        jxz1=ixz1                                                       2d26s20
        jyz1=iyz1                                                       2d26s20
        jxy1=ixy1                                                       2d26s20
        do id=0,nds-1                                                    2d24s20
         idp=id+ibsstor(ipack2(4))+nbaslarge*2                          2d25s20
         idpp=idp+nbassmall                                             2d25s20
         do ic=0,ncl-1                                                   2d24s20
          icp=ic+iblstor(ipack2(3))                                      2d24s20
          icpp=icp+nbaslarge                                             2d24s20
          jfockaar=icp+1+nbtot*idp                                         2d24s20
          jfockaai=jfockaar+nb2                                            2d25s20
          jfockabr=icp+1+nbtot*idpp                                         2d24s20
          jfockabi=jfockabr+nb2                                            2d25s20
          jfockbar=icpp+1+nbtot*idp                                     2d25s20
          jfockbai=jfockbar+nb2                                            2d24s20
          jfockbbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbbi=jfockbbr+nb2                                            2d24s20
          do ib=0,nbl-1                                                  2d24s20
           ibp=ib+iblstor(ipack2(2))                                    2d26s20
           ibpp=ibp+nbaslarge                                           2d26s20
           do ia=0,nas-1                                                 2d24s20
            iap=ia+ibsstor(ipack2(1))+nbaslarge*2                       2d26s20
            iapp=iap+nbassmall                                          2d26s20
            jdenaar=iap+1+nbtot*ibp                                     2d25s20
            jdenaai=jdenaar+nb2                                         2d25s20
            jdenabr=iap+1+nbtot*ibpp                                    2d25s20
            jdenabi=jdenabr+nb2                                         2d25s20
            jdenbar=iapp+1+nbtot*ibp                                    2d25s20
            jdenbai=jdenbar+nb2                                         2d25s20
            jdenbbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbbi=jdenbbr+nb2                                            2d24s20
            termr=bc(jzz1)*(den(jdenaar)-den(jdenbbr))                  2d26s20
     $           +bc(jxz1)*(den(jdenabr)+den(jdenbar))                  2d26s20
     $           +bc(jyz1)*(-den(jdenbai)+den(jdenabi))                 2d26s20
            termi=bc(jzz1)*(den(jdenaai)-den(jdenbbi))                  2d26s20
     $           +bc(jxz1)*(den(jdenabi)+den(jdenbai))                  2d26s20
     $           +bc(jyz1)*(den(jdenbar)-den(jdenabr))                  2d26s20
            fock(jfockaar)=fock(jfockaar)+termr                         2d25s20
            fock(jfockaai)=fock(jfockaai)+termi                         2d25s20
            fock(jfockbbr)=fock(jfockbbr)-termr                         2d26s20
            fock(jfockbbi)=fock(jfockbbi)-termi                         2d26s20
            ddr=den(jdenaar)-den(jdenbbr)                               2d26s20
            ddi=den(jdenaai)-den(jdenbbi)                               2d26s20
            termr1=bc(jxz1)*ddr+bc(jxx1)*(den(jdenbar)+den(jdenabr))
            termr2=-bc(jyz1)*ddi                                        2d26s20
            termi1=bc(jxz1)*ddi+bc(jxx1)*(den(jdenbai)+den(jdenabi))
            termi2=bc(jyz1)*ddr
            fock(jfockabr)=fock(jfockabr)+termr1-termr2
     $           +2d0*bc(jxy1)*den(jdenabi)
            fock(jfockabi)=fock(jfockabi)+termi1-termi2
     $           -2d0*bc(jxy1)*den(jdenabr)
            fock(jfockbar)=fock(jfockbar)+termr1+termr2
     $           -2d0*bc(jxy1)*den(jdenbai)
            fock(jfockbai)=fock(jfockbai)+termi1+termi2
     $           +2d0*bc(jxy1)*den(jdenbar)
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            termr=bc(jxz1)*(den(kdenbar)+den(kdenabr))
     $           -bc(jyz1)*(den(kdenbai)-den(kdenabi))
            termi=bc(jxz1)*(den(kdenbai)+den(kdenabi))                  2d26s20
     $           +bc(jyz1)*(den(kdenbar)-den(kdenabr))                  2d26s20
            fock(kfockaar)=fock(kfockaar)-termr-den(kdenaar)*bc(jzz1)   2d26s20
     $           -den(kdenbbr)*bc(jxx1)                                 2d26s20
            fock(kfockaai)=fock(kfockaai)-termi-den(kdenaai)*bc(jzz1)   2d26s20
     $           -den(kdenbbi)*bc(jxx1)                                 2d26s20
            fock(kfockbbr)=fock(kfockbbr)+termr-den(kdenaar)*bc(jxx1)   2d26s20
     $           -den(kdenbbr)*bc(jzz1)                                 2d26s20
            fock(kfockbbi)=fock(kfockbbi)+termi-den(kdenaai)*bc(jxx1)   2d26s20
     $           -den(kdenbbi)*bc(jzz1)                                 2d26s20
            ddr=den(kdenaar)-den(kdenbbr)                               2d26s20
            ddi=den(kdenaai)-den(kdenbbi)                               2d26s20
            termr1=bc(jxz1)*ddr                                         2d26s20
            termr2=-bc(jyz1)*ddi                                        2d26s20
            termi1=bc(jxz1)*ddi                                         2d26s20
            termi2=bc(jyz1)*ddr
            fock(kfockabr)=fock(kfockabr)-termr1+termr2
     $           -bc(jxx1)*den(kdenabr)+bc(jzz1)*den(kdenbar)
     $           -2d0*bc(jxy1)*den(kdenabi)                             2d26s20
            fock(kfockabi)=fock(kfockabi)-termi1+termi2
     $           -bc(jxx1)*den(kdenabi)+bc(jzz1)*den(kdenbai)
     $           +2d0*bc(jxy1)*den(kdenabr)                             2d26s20
            fock(kfockbar)=fock(kfockbar)-termr1-termr2
     $           -bc(jxx1)*den(kdenbar)+bc(jzz1)*den(kdenabr)
     $           +2d0*bc(jxy1)*den(kdenbai)                             2d26s20
            fock(kfockbai)=fock(kfockbai)-termi1-termi2
     $           -bc(jxx1)*den(kdenbai)+bc(jzz1)*den(kdenabi)
     $           -2d0*bc(jxy1)*den(kdenbar)                             2d26s20
            jzz1=jzz1+1                                                    2d24s20
            jxx1=jxx1+1                                                    2d24s20
            jxz1=jxz1+1                                                    2d24s20
            jyz1=jyz1+1                                                    2d24s20
            jxy1=jxy1+1                                                    2d24s20
           end do
          end do
         end do
        end do
        ibcoff=inorm                                                    2d26s20
        nla=2*la8+1                                                     2d25s20
        nlb=2*ibc(jb)+1                                                 2d26s20
        nlc=2*lc8+1                                                     2d26s20
        nld=2*ibc(jd)+1                                                 2d26s20
        xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $       *dfloat(irtyp(5))                                          2d26s20
        xnorm=-xnorm*0.5d0*scaleb                                              2d26s20
        nabcd=nal*nbs*ncl*nds                                           2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,inorm,0,0,0,0,0,0,0,0,0,0,0,0,                    2d26s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=inorm+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,izz1,0,0,1,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=izz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,izz2,0,0,0,0,0,1,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(izz1+ii)=bc(izz1+ii)+bc(izz2+ii)+bc(inorm+ii)               2d26s20
        end do                                                          2d26s20
        ibcoff=izz2                                                     2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixx1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         2,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixx1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixx2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         2,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixx2+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyy1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=iyy1+nabcd
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyy2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixx1+ii)=bc(ixx1+ii)+bc(ixx2+ii)+bc(iyy1+ii)+bc(iyy2+ii)    2d26s20
     $        +2d0*bc(inorm+ii)                                         2d26s20
        end do
        ibcoff=ixx2
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixz1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixz2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixz1+ii)=bc(ixz1+ii)+bc(ixz2+ii)                               2d26s20
        end do                                                          2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyz1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=iyz1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,iyz2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $         4,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(iyz1+ii)=bc(iyz1+ii)+bc(iyz2+ii)                               2d26s20
        end do                                                          2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixy1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        ibcoff=ixy1+nabcd                                               2d26s20
        call derid(                                                      2d24s20
     $      la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ixy2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $         3,.false.,xnorm,bc,ibc)                                  11d17s22
        do ii=0,nabcd-1                                                  2d26s20
         bc(ixy1+ii)=bc(ixy1+ii)+bc(ixy2+ii)                               2d26s20
        end do                                                          2d26s20
        jzz1=izz1                                                       2d26s20
        jxx1=ixx1                                                       2d26s20
        jxz1=ixz1                                                       2d26s20
        jyz1=iyz1                                                       2d26s20
        jxy1=ixy1                                                       2d26s20
        do id=0,ndl-1                                                    2d24s20
         idp=id+iblstor(ipack2(4))                                      2d26s20
         idpp=idp+nbaslarge                                             2d26s20
         do ic=0,ncs-1                                                   2d24s20
          icp=ic+ibsstor(ipack2(3))+nbaslarge*2                         2d26s20
          icpp=icp+nbassmall                                            2d26s20
          jfockaar=icp+1+nbtot*idp                                         2d24s20
          jfockaai=jfockaar+nb2                                            2d25s20
          jfockabr=icp+1+nbtot*idpp                                         2d24s20
          jfockabi=jfockabr+nb2                                            2d25s20
          jfockbar=icpp+1+nbtot*idp                                     2d25s20
          jfockbai=jfockbar+nb2                                            2d24s20
          jfockbbr=icpp+1+nbtot*idpp                                      2d24s20
          jfockbbi=jfockbbr+nb2                                            2d24s20
          do ib=0,nbl-1                                                  2d24s20
           ibp=ib+iblstor(ipack2(2))                                    2d26s20
           ibpp=ibp+nbaslarge                                           2d26s20
           do ia=0,nas-1                                                 2d24s20
            iap=ia+ibsstor(ipack2(1))+nbaslarge*2                       2d26s20
            iapp=iap+nbassmall                                          2d26s20
            jdenaar=iap+1+nbtot*ibp                                     2d25s20
            jdenaai=jdenaar+nb2                                         2d25s20
            jdenabr=iap+1+nbtot*ibpp                                    2d25s20
            jdenabi=jdenabr+nb2                                         2d25s20
            jdenbar=iapp+1+nbtot*ibp                                    2d25s20
            jdenbai=jdenbar+nb2                                         2d25s20
            jdenbbr=iapp+1+nbtot*ibpp                                     2d25s20
            jdenbbi=jdenbbr+nb2                                            2d24s20
            termr=bc(jzz1)*(den(jdenaar)-den(jdenbbr))                  2d26s20
     $           +bc(jxz1)*(den(jdenabr)+den(jdenbar))                  2d26s20
     $           +bc(jyz1)*(-den(jdenbai)+den(jdenabi))                 2d26s20
            termi=bc(jzz1)*(den(jdenaai)-den(jdenbbi))                  2d26s20
     $           +bc(jxz1)*(den(jdenabi)+den(jdenbai))                  2d26s20
     $           +bc(jyz1)*(den(jdenbar)-den(jdenabr))                  2d26s20
            fock(jfockaar)=fock(jfockaar)+termr                         2d25s20
            fock(jfockaai)=fock(jfockaai)+termi                         2d25s20
            fock(jfockbbr)=fock(jfockbbr)-termr                         2d26s20
            fock(jfockbbi)=fock(jfockbbi)-termi                         2d26s20
            ddr=den(jdenaar)-den(jdenbbr)                               2d26s20
            ddi=den(jdenaai)-den(jdenbbi)                               2d26s20
            termr1=bc(jxz1)*ddr+bc(jxx1)*(den(jdenbar)+den(jdenabr))
            termr2=-bc(jyz1)*ddi                                        2d26s20
            termi1=bc(jxz1)*ddi+bc(jxx1)*(den(jdenbai)+den(jdenabi))
            termi2=bc(jyz1)*ddr
            fock(jfockabr)=fock(jfockabr)+termr1-termr2
     $           +2d0*bc(jxy1)*den(jdenabi)
            fock(jfockabi)=fock(jfockabi)+termi1-termi2
     $           -2d0*bc(jxy1)*den(jdenabr)
            fock(jfockbar)=fock(jfockbar)+termr1+termr2
     $           -2d0*bc(jxy1)*den(jdenbai)
            fock(jfockbai)=fock(jfockbai)+termi1+termi2
     $           +2d0*bc(jxy1)*den(jdenbar)
            kdenaar=iap+1+nbtot*idp                                        2d24s20
            kdenaai=kdenaar+nb2                                              2d24s20
            kdenbar=iapp+1+nbtot*idp                                     2d25s20
            kdenbai=kdenbar+nb2                                          2d25s20
            kdenabr=iap+1+nbtot*idpp                                     2d25s20
            kdenabi=kdenabr+nb2                                          2d25s20
            kdenbbr=iapp+1+nbtot*idpp                                     2d24s20
            kdenbbi=kdenbbr+nb2                                            2d24s20
            kfockaar=icp+1+nbtot*ibp                                       2d24s20
            kfockaai=kfockaar+nb2                                            2d24s20
            kfockabr=icp+1+nbtot*ibpp                                    2d25s20
            kfockabi=kfockabr+nb2                                        2d25s20
            kfockbar=icpp+1+nbtot*ibp                                    2d25s20
            kfockbai=kfockbar+nb2                                        2d25s20
            kfockbbr=icpp+1+nbtot*ibpp                                    2d24s20
            kfockbbi=kfockbbr+nb2                                          2d24s20
            termr=bc(jxz1)*(den(kdenbar)+den(kdenabr))
     $           -bc(jyz1)*(den(kdenbai)-den(kdenabi))
            termi=bc(jxz1)*(den(kdenbai)+den(kdenabi))                  2d26s20
     $           +bc(jyz1)*(den(kdenbar)-den(kdenabr))                  2d26s20
            fock(kfockaar)=fock(kfockaar)-termr-den(kdenaar)*bc(jzz1)   2d26s20
     $           -den(kdenbbr)*bc(jxx1)                                 2d26s20
            fock(kfockaai)=fock(kfockaai)-termi-den(kdenaai)*bc(jzz1)   2d26s20
     $           -den(kdenbbi)*bc(jxx1)                                 2d26s20
            fock(kfockbbr)=fock(kfockbbr)+termr-den(kdenaar)*bc(jxx1)   2d26s20
     $           -den(kdenbbr)*bc(jzz1)                                 2d26s20
            fock(kfockbbi)=fock(kfockbbi)+termi-den(kdenaai)*bc(jxx1)   2d26s20
     $           -den(kdenbbi)*bc(jzz1)                                 2d26s20
            ddr=den(kdenaar)-den(kdenbbr)                               2d26s20
            ddi=den(kdenaai)-den(kdenbbi)                               2d26s20
            termr1=bc(jxz1)*ddr                                         2d26s20
            termr2=-bc(jyz1)*ddi                                        2d26s20
            termi1=bc(jxz1)*ddi                                         2d26s20
            termi2=bc(jyz1)*ddr
            fock(kfockabr)=fock(kfockabr)-termr1+termr2
     $           -bc(jxx1)*den(kdenabr)+bc(jzz1)*den(kdenbar)
     $           -2d0*bc(jxy1)*den(kdenabi)                             2d26s20
            fock(kfockabi)=fock(kfockabi)-termi1+termi2
     $           -bc(jxx1)*den(kdenabi)+bc(jzz1)*den(kdenbai)
     $           +2d0*bc(jxy1)*den(kdenabr)                             2d26s20
            fock(kfockbar)=fock(kfockbar)-termr1-termr2
     $           -bc(jxx1)*den(kdenbar)+bc(jzz1)*den(kdenabr)
     $           +2d0*bc(jxy1)*den(kdenbai)                             2d26s20
            fock(kfockbai)=fock(kfockbai)-termi1-termi2
     $           -bc(jxx1)*den(kdenbai)+bc(jzz1)*den(kdenabi)
     $           -2d0*bc(jxy1)*den(kdenbar)                             2d26s20
            jzz1=jzz1+1                                                    2d24s20
            jxx1=jxx1+1                                                    2d24s20
            jxz1=jxz1+1                                                    2d24s20
            jyz1=jyz1+1                                                    2d24s20
            jxy1=jxy1+1                                                    2d24s20
           end do
          end do
         end do
        end do
        ibcoff=inorm                                                    2d26s20
       end if                                                           2d26s20
      end do                                                            2d19s10
      call dws_gsumf(fock,nbb)                                           5d3s10
      ibcoff=ibcoffo                                                    2d19s10
      return
      end                                                               2d19s10
