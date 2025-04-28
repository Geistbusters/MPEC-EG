c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine paraeri4ct(natom,ngaus,ibdat,nbasis,oooo,nox,vec,      2d27s20
     $                iblstor,ibsstor,idwsdeb,ascale,nbaslarge,         2d28s20
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
      dimension oooo(*),vec(nbtot,nbtot,2),iblstor(*),ibsstor(*),       2d28s20
     $     irtyp(5)                                                     2d28s20
      scaleb=1d0
      nbaslarge2=nbaslarge*2                                            2d29s20
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
      nnt=ngaus**2                                                      2d27s20
      ipair=ibcoff                                                      2d19s10
      ibcoff=ipair+nnt                                                  2d24s20
      i12=ibcoff                                                        2d19s10
      ibcoff=i12+nnt*2                                                   2d19s10
      call enough('paraeri4ct.  1',bc,ibc)
      j12=i12                                                           2d19s10
      jbdat=ibdat-1                                                     2d19s10
      ngaus7=ngaus*7                                                    5d3s10
      ii12=0
      do i1=1,ngaus                                                     2d19s10
       ipack2(1)=i1
       do i2=1,ngaus
        i12sz=((ibc(jbdat+i1)+ibc(jbdat+i2)+2)**2)                      2d24s20
        ipack2(2)=i2                                                    2d24s20
        ibc(j12)=-i12sz                                                 2d27s20
        ibc(j12+nnt)=ipack8                                             2d27s20
        j12=j12+1
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
      call second(time2)
      telap=time2-time1
      loopit=0                                                          5d25s18
      ihalfr=ibcoff                                                     2d28s20
      nwdsh=nbtot*nbtot*nox*nox                                         2d28s20
      ihalfi=ihalfr+nwdsh                                               2d28s20
      ibcoff=ihalfi+nwdsh                                               2d28s20
      call enough('paraeri4ct.  2',bc,ibc)
      do i=ihalfr,ibcoff-1                                              2d28s20
       bc(i)=0d0                                                        2d28s20
      end do                                                            2d28s20
      do i=1+mynowprog,nnt,mynprocg                                     2d24s20
       loopit=loopit+1                                                  5d25s18
       jpair=ipair+i-1                                                  2d24s20
       ipack8=ibc(jpair)                                                2d24s20
       ishella=ipack2(1)                                                2d28s20
       ishellb=ipack2(2)                                                2d28s20
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
       nal=ibc(ilsz+ibc(ja))                                            2d24s20
       nbl=ibc(ilsz+ibc(jb))                                            2d24s20
       nas=ibc(ilsz+ibc(ja)+1)                                          2d24s20
       nbs=ibc(ilsz+ibc(jb)+1)                                          2d24s20
       la8=ibc(ja)+1                                                    2d24s20
       lb8=ibc(jb)+1                                                    2d24s20
       itmpe=ibcoff                                                     2d27s20
       naa=2*(nal+nas)                                                  2d27s20
       nbb=2*(nbl+nbs)                                                  2d27s20
       nwds=nbtot*nbtot*naa*nbb                                         2d27s20
       ibcoff=itmpe+nwds                                                2d27s20
       if(irtyp(5).ne.0)then                                            2d27s20
        itmpei=ibcoff                                                   2d27s20
        ibcoff=itmpei+nwds                                              2d27s20
       end if                                                           2d27s20
       call enough('paraeri4ct.  3',bc,ibc)
       do ii=itmpe,ibcoff-1                                              2d27s20
        bc(ii)=0d0                                                       2d27s20
       end do                                                           2d27s20
       do icd=0,nnt-1                                                   2d27s20
        ipack8=ibc(ipair+icd)                                           2d27s20
        ishellc=ipack2(1)                                               3d28s20
        ishelld=ipack2(2)                                                2d28s20
        jc=jbdat+ipack2(1)
        isnormc=isnorm+ipack2(1)-1
        jc2=jc+ngaus
        jc3=jc2+ngaus
        jc4=jc3+ngaus
        jc5=jc4+ngaus
        jc6=jc5+ngaus
        jc7=jc6+ngaus
        jc8=jc7+ngaus
        jd=jbdat+ipack2(2)
        isnormd=isnorm+ipack2(2)-1
        jd2=jd+ngaus
        jd3=jd2+ngaus
        jd4=jd3+ngaus
        jd5=jd4+ngaus
        jd6=jd5+ngaus
        jd7=jd6+ngaus
        jd8=jd7+ngaus
        ncl=ibc(ilsz+ibc(jc))                                            2d24s20
        ndl=ibc(ilsz+ibc(jd))                                            2d24s20
        ncs=ibc(ilsz+ibc(jc)+1)                                          2d24s20
        nds=ibc(ilsz+ibc(jd)+1)                                          2d24s20
        lc8=ibc(jc)+1                                                    2d24s20
        ld8=ibc(jd)+1                                                    2d24s20
c
c     shadow to get symmetry right
c
        kc=jd
        kd=jc
        ksnormc=isnormd
        ksnormd=isnormc
        kc2=jd2
        kd2=jc2
        kc3=jd3
        kd3=jc3
        kc4=jd4
        kd4=jc4
        kc5=jd5
        kd5=jc5
        kc6=jd6
        kd6=jc6
        kc7=jd7
        kd7=jc7
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
         jb1=ib1                                                          2d24s20
         do id=0,ndl-1                                                    2d24s20
          idp=id+iblstor(ishelld)                                       3d28s20
          idpp=idp+nbaslarge                                              2d24s20
          do ic=0,ncl-1                                                   2d24s20
           icp=ic+iblstor(ishellc)                                      3d28s20
           icpp=icp+nbaslarge                                             2d24s20
           do ib=0,nbl-1                                                  2d24s20
            ibp=ib+nbl                                                   2d27s20
            do ia=0,nal-1                                                 2d24s20
             iap=ia+nal                                                  2d27s20
             iad1=itmpe+ia+naa*(ib+nbb*(icp+nbtot*idp))                  2d27s20
             iad2=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))                2d27s20
             iad3=itmpe+ia+naa*(ib+nbb*(icpp+nbtot*idpp))                2d27s20
             iad4=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))              2d27s20
             bc(iad1)=bc(jb1)                                            2d27s20
             bc(iad2)=bc(jb1)                                            2d27s20
             bc(iad3)=bc(jb1)                                            2d27s20
             bc(iad4)=bc(jb1)                                            2d27s20
             jb1=jb1+1                                                    2d24s20
            end do
           end do
          end do
         end do
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
          idp=id+iblstor(ishelld)                                       3d28s20
          idpp=idp+nbaslarge                                              2d24s20
          do ic=0,ncl-1                                                   2d24s20
           icp=ic+iblstor(ishellc)                                      3d28s20
           icpp=icp+nbaslarge                                             2d24s20
           do ib=0,nbs-1                                                  2d24s20
            lbp=ib+ibsstor(ishellb)+nbaslarge2                          3d27s20
            lbpp=lbp+nbassmall                                          3d27s20
            ibp=ib+2*nbl                                                2d27s20
            ibpp=ibp+nbs                                                2d27s20
            do ia=0,nas-1                                                 2d24s20
             iap=ia+2*nal                                                2d27s20
             iapp=iap+nas                                                2d27s20
             lap=ia+ibsstor(ishella)+nbaslarge2                         3d27s20
             lapp=lap+nbassmall                                         3d27s20
             iad1=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))               2d27s20
             iad2=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))               2d27s20
             iad3=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))             2d27s20
             iad4=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))             2d27s20
             bc(iad1)=bc(jb1)                                           2d27s20
             bc(iad2)=bc(jb1)                                           2d27s20
             bc(iad3)=bc(jb1)                                           2d27s20
             bc(iad4)=bc(jb1)                                           2d27s20
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
          idp=id+ibsstor(ishelld)+nbaslarge2                            3d28s20
          idpp=idp+nbassmall                                             2d25s20
          do ic=0,ncs-1                                                   2d24s20
           icp=ic+ibsstor(ishellc)+nbaslarge2                           3d28s20
           icpp=icp+nbassmall                                            2d25s20
           do ib=0,nbl-1                                                  2d24s20
            ibp=ib+nbl                                                  2d27s20
            lbp=ib+iblstor(ishellb)                                     3d27s20
            llbp=lbp+nbaslarge                                          3d27s20
            do ia=0,nal-1                                                 2d24s20
             iap=ia+nal                                                 2d27s20
             lap=ia+iblstor(ishella)                                    3d27s20
             lapp=lap+nbaslarge                                         3d27s20
             iad1=itmpe+ia+naa*(ib+nbb*(icp+nbtot*idp))                 2d27s20
             iad2=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))                 2d27s20
             iad3=itmpe+ia+naa*(ib+nbb*(icpp+nbtot*idpp))                 2d27s20
             iad4=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))                 2d27s20
             bc(iad1)=bc(jb1)
             bc(iad2)=bc(jb1)
             bc(iad3)=bc(jb1)
             bc(iad4)=bc(jb1)
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
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d24s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        1,.false.,xnorm,bc,ibc)                                   11d17s22
         jb1=ib1                                                          2d24s20
         do id=0,nds-1                                                    2d24s20
          idp=id+ibsstor(ipack2(2))+nbaslarge2                          2d29s20
          idpp=idp+nbassmall                                             2d25s20
          do ic=0,ncs-1                                                   2d24s20
           icp=ic+ibsstor(ipack2(1))+nbaslarge2                         2d29s20
           icpp=icp+nbassmall                                            2d25s20
           do ib=0,nbs-1                                                  2d24s20
            ibp=ib+nbl*2                                                2d27s20
            ibpp=ibp+nbs                                                2d27s20
            do ia=0,nas-1                                                 2d24s20
             iap=ia+nal*2                                               2d27s20
             iapp=iap+nas                                               2d27s20
             iad1=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))               2d27s20
             iad2=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))               2d27s20
             iad3=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))               2d27s20
             iad4=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))           2d27s20
             bc(iad1)=bc(jb1)                                           2d27s20
             bc(iad2)=bc(jb1)                                           2d27s20
             bc(iad3)=bc(jb1)                                           2d27s20
             bc(iad4)=bc(jb1)                                           2d27s20
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
     $        *dfloat(irtyp(4))                                          2d25s20
         xnorm=-xnorm*scaleb
         if(irtyp(5).eq.0)xnorm=xnorm*0.5d0                             3d28s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        1,.false.,xnorm,bc,ibc)                                   11d17s22
         jb1=ib1                                                          2d24s20
         do id=0,nds-1                                                    2d24s20
          idp=id+ibsstor(ipack2(2))+nbaslarge2                           2d25s20
          idpp=idp+nbassmall                                             2d25s20
          do ic=0,ncl-1                                                   2d24s20
           icp=ic+iblstor(ipack2(1))                                      2d24s20
           icpp=icp+nbaslarge                                             2d24s20
           do ib=0,nbs-1                                                  2d24s20
            ibp=ib+2*nbl
            ibpp=ibp+nbs                                                2d27s20
            do ia=0,nal-1                                                 2d24s20
             iap=ia                                                     2d27s20
             iapp=iap+nal                                                 2d27s20
             iad11=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))                2d27s20
             iad22=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))                2d27s20
             iad33=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))                2d27s20
             iad44=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))                2d27s20
             iad23=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idpp))                2d27s20
             iad32=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idp))                2d27s20
             bc(iad11)=bc(jb1)                                          2d27s20
             bc(iad22)=-bc(jb1)                                          2d27s20
             bc(iad33)=-bc(jb1)                                         2d28s20
             bc(iad44)=bc(jb1)                                          2d28s20
             bc(iad23)=2d0*bc(jb1)                                          2d27s20
             bc(iad32)=2d0*bc(jb1)                                          2d27s20
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
         if(irtyp(5).eq.0)xnorm=xnorm*0.5d0                             3d28s20
         call derid(                                                      2d24s20
     $      ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $      lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $      lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $      ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $      dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
         jb1=ib1                                                          2d24s20
         do id=0,ndl-1                                                    2d24s20
          idp=id+iblstor(ipack2(2))                                      2d25s20
          idpp=idp+nbaslarge                                             2d25s20
          do ic=0,ncs-1                                                   2d24s20
           icp=ic+ibsstor(ipack2(1))+nbaslarge2                          2d25s20
           icpp=icp+nbassmall                                            2d25s20
           do ib=0,nbs-1                                                  2d24s20
            ibp=ib+2*nbl                                                2d27s20
            ibpp=ibp+nbs                                                2d27s20
            do ia=0,nal-1                                                 2d24s20
             iap=ia
             iapp=iap+nal                                                 2d27s20
             iad11=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))                2d27s20
             iad22=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))                2d27s20
             iad33=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))                2d27s20
             iad44=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))                2d27s20
             iad23=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idpp))                2d27s20
             iad32=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idp))                2d27s20
             bc(iad11)=bc(jb1)                                          2d27s20
             bc(iad22)=-bc(jb1)                                          2d27s20
             bc(iad33)=-bc(jb1)                                         2d28s20
             bc(iad44)=bc(jb1)                                          2d28s20
             bc(iad23)=2d0*bc(jb1)                                          2d27s20
             bc(iad32)=2d0*bc(jb1)                                          2d27s20
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
     $        *dfloat(irtyp(4))                                          2d25s20
         xnorm=-xnorm*scaleb
         if(irtyp(5).eq.0)xnorm=xnorm*0.5d0                             3d28s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
         jb1=ib1                                                          2d24s20
         do id=0,nds-1                                                    2d24s20
          idp=id+ibsstor(ipack2(2))+nbaslarge2                          2d29s20
          idpp=idp+nbassmall                                             2d25s20
          do ic=0,ncl-1                                                   2d24s20
           icp=ic+iblstor(ipack2(1))                                      2d24s20
           icpp=icp+nbaslarge                                             2d24s20
           do ib=0,nbl-1                                                  2d24s20
            ibp=ib
            ibpp=ibp+nbl                                                  2d27s20
            do ia=0,nas-1                                                 2d24s20
             iap=ia+2*nal                                               2d27s20
             iapp=iap+nas                                               2d27s20
             iad11=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))                2d27s20
             iad22=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))                2d27s20
             iad33=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))                2d27s20
             iad44=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))                2d27s20
             iad23=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idpp))                2d27s20
             iad32=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idp))                2d27s20
             bc(iad11)=bc(jb1)                                          2d27s20
             bc(iad22)=-bc(jb1)                                          2d27s20
             bc(iad33)=-bc(jb1)                                         2d28s20
             bc(iad44)=bc(jb1)                                          2d28s20
             bc(iad23)=2d0*bc(jb1)                                          2d27s20
             bc(iad32)=2d0*bc(jb1)                                          2d27s20
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
     $        *dfloat(irtyp(4))                                          2d25s20
         xnorm=-xnorm*scaleb
         if(irtyp(5).eq.0)xnorm=xnorm*0.5d0                             3d28s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d25s20
     $        dum,1,idum,ib1,0,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $         1,.false.,xnorm,bc,ibc)                                  11d17s22
         jb1=ib1                                                          2d24s20
         do id=0,ndl-1                                                    2d24s20
          idp=id+iblstor(ipack2(2))                                      2d25s20
          idpp=idp+nbaslarge                                             2d25s20
          do ic=0,ncs-1                                                   2d24s20
           icp=ic+ibsstor(ipack2(1))+nbaslarge2                         2d29s20
           icpp=icp+nbassmall                                            2d25s20
           do ib=0,nbl-1                                                  2d24s20
            ibp=ib
            ibpp=ibp+nbl
            do ia=0,nas-1                                                 2d24s20
             iap=ia+2*nal                                               2d27s20
             iapp=iap+nas                                               2d27s20
             iad11=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))                2d27s20
             iad22=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))                2d27s20
             iad33=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))                2d27s20
             iad44=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))                2d27s20
             iad23=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idpp))                2d27s20
             iad32=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idp))                2d27s20
             bc(iad11)=bc(jb1)                                          2d27s20
             bc(iad22)=-bc(jb1)                                          2d27s20
             bc(iad33)=-bc(jb1)                                         2d28s20
             bc(iad44)=bc(jb1)                                          2d28s20
             bc(iad23)=2d0*bc(jb1)                                          2d27s20
             bc(iad32)=2d0*bc(jb1)                                          2d27s20
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
         fm=1d0
         fo=1d0
         nla=2*ibc(ja)+1                                                     2d25s20
         nlb=2*lb8+1                                                     2d25s20
         nlc=2*ibc(jc)+1                                                  2d24s20
         nld=2*ld8+1                                                     2d25s20
         xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $        *dfloat(irtyp(5))                                          2d26s20
         xnorm=-xnorm*0.5d0*scaleb                                              2d26s20
         nabcd=nal*nbs*ncl*nds                                           2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,izz1,0,0,1,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=izz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,izz2,0,0,0,0,0,1,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(izz1+ii)=fo*(bc(izz1+ii)+bc(izz2+ii))
         end do                                                          2d26s20
         ibcoff=izz2                                                     2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixx1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        2,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixx1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixx2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        2,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixx2+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyy1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=iyy1+nabcd
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyy2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          xsum=bc(ixx1+ii)+bc(ixx2+ii)                                  2d29s20
          ysum=bc(iyy1+ii)+bc(iyy2+ii)                                  2d29s20
          bc(ixx1+ii)=fo*(xsum+ysum)
          bc(ixx2+ii)=fo*(xsum-ysum)                                         2d29s20
         end do
         ibcoff=iyy1                                                    2d29s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixz1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixz2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(ixz1+ii)=fo*(bc(ixz1+ii)+bc(ixz2+ii))                               2d26s20
         end do                                                          2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyz1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=iyz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyz2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(iyz1+ii)=fo*(bc(iyz1+ii)+bc(iyz2+ii))                               2d26s20
         end do                                                          2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixy1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixy1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixy2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(ixy1+ii)=fo*(bc(ixy1+ii)+bc(ixy2+ii))                               2d26s20
         end do                                                          2d26s20
         jzz1=izz1                                                       2d26s20
         jxx1=ixx1                                                       2d26s20
         jxx2=ixx2                                                      3d1s20
         jxz1=ixz1                                                       2d26s20
         jyz1=iyz1                                                       2d26s20
         jxy1=ixy1                                                       2d26s20
         do id=0,nds-1                                                    2d24s20
          idp=id+ibsstor(ipack2(2))+nbaslarge2                          2d29s20
          idpp=idp+nbassmall                                             2d25s20
          do ic=0,ncl-1                                                   2d24s20
           icp=ic+iblstor(ipack2(1))                                      2d24s20
           icpp=icp+nbaslarge                                             2d24s20
           do ib=0,nbs-1                                                  2d24s20
            ibp=ib+2*nbl                                                2d27s20
            ibpp=ibp+nbs                                                2d27s20
            do ia=0,nal-1                                                 2d24s20
             iap=ia
             iapp=iap+nal                                                 2d27s20
             iad11=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))              2d27s20 aaaa     4
             iad22=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))              2d27s20 bbaa   2 & 2 sym
             iad33=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))              2d27s20 aabb   2 & 2 sym
             iad44=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))              2d27s20 bbbb 4
             bc(iad11)=bc(iad11)+bc(jzz1)
             bc(iad22)=bc(iad22)-bc(jzz1)
             bc(iad33)=bc(iad33)-bc(jzz1)
             bc(iad44)=bc(iad44)+bc(jzz1)
             iad41=itmpe+iapp+naa*(ibp+nbb*(icpp+nbtot*idp))            2d27s20 baba
             iad32=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idp))            2d27s20 abba  2 & 2, asym
             iad23=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idpp))            2d27s20 baab
             iad14=itmpe+iap+naa*(ibpp+nbb*(icp+nbtot*idpp))            2d27s20 abab
             bc(iad41)=bc(jxx2)                                         3d1s20
             bc(iad32)=bc(iad32)+bc(jxx1)                                         2d27s20
             bc(iad23)=bc(iad23)+bc(jxx1)                                         2d27s20
             bc(iad14)=bc(jxx2)                                         3d1s20
             iad41=iad41+nwds                                           2d27s20
             iad14=iad14+nwds                                           2d27s20
             bc(iad41)=2d0*bc(jxy1)                                     2d27s20
             bc(iad14)=-2d0*bc(jxy1)                                    2d27s20
             iad21=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idp))             2d27s20 baaa 1b
             iad31=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idp))             2d27s20 aaba
             iad12=itmpe+iap+naa*(ibpp+nbb*(icp+nbtot*idp))             2d27s20 abaa
             iad13=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idpp))             2d27s20 aaab
             bc(iad21)=bc(jxz1)                                         2d27s20
             bc(iad31)=bc(jxz1)                                         2d27s20
             bc(iad12)=bc(jxz1)                                         2d27s20
             bc(iad13)=bc(jxz1)                                         2d27s20
             iad21=iad21+nwds                                           2d27s20
             iad31=iad31+nwds                                           2d27s20
             iad12=iad12+nwds                                           2d27s20
             iad13=iad13+nwds                                           2d27s20
             bc(iad21)=bc(jyz1)                                         2d27s20
             bc(iad31)=bc(jyz1)                                         2d27s20
             bc(iad12)=-bc(jyz1)                                         2d27s20
             bc(iad13)=-bc(jyz1)                                         2d27s20
             iad42=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idp))           2d27s20 bbba     1a
             iad43=itmpe+iapp+naa*(ibp+nbb*(icpp+nbtot*idpp))           2d27s20 babb
             iad24=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idpp))           2d27s20 bbab
             iad34=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idpp))           2d27s20 abbb
             bc(iad42)=-bc(jxz1)                                        2d27s20
             bc(iad43)=-bc(jxz1)                                        2d27s20
             bc(iad24)=-bc(jxz1)                                        2d27s20
             bc(iad34)=-bc(jxz1)                                        2d27s20
             iad42=iad42+nwds                                           2d27s20
             iad43=iad43+nwds                                           2d27s20
             iad24=iad24+nwds                                           2d27s20
             iad34=iad34+nwds                                           2d27s20
             bc(iad42)=-bc(jyz1)                                        2d27s20
             bc(iad43)=-bc(jyz1)                                        2d27s20
             bc(iad24)=bc(jyz1)                                         2d27s20
             bc(iad34)=bc(jyz1)                                         2d27s20
             jzz1=jzz1+1                                                    2d24s20
             jxx1=jxx1+1                                                    2d24s20
             jxx2=jxx2+1                                                3d1s20
             jxz1=jxz1+1                                                    2d24s20
             jyz1=jyz1+1                                                    2d24s20
             jxy1=jxy1+1                                                    2d24s20
            end do
           end do
          end do
         end do
         ibcoff=izz1                                                    3d1s20
         nla=2*ibc(ja)+1                                                     2d25s20
         nlb=2*lb8+1                                                     2d25s20
         nlc=2*lc8+1                                                     2d26s20
         nld=2*ibc(jd)+1                                                 2d26s20
         xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $        *dfloat(irtyp(5))                                          2d26s20
         xnorm=-xnorm*0.5d0*scaleb                                              2d26s20
         nabcd=nal*nbs*ncs*ndl                                          3d1s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,izz1,0,0,1,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=izz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,izz2,0,0,0,0,0,1,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(izz1+ii)=fo*(bc(izz1+ii)+bc(izz2+ii))
         end do                                                          2d26s20
         ibcoff=izz2                                                     2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixx1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        2,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixx1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixx2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        2,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixx2+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyy1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=iyy1+nabcd
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyy2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          xsum=bc(ixx1+ii)+bc(ixx2+ii)                                  2d29s20
          ysum=bc(iyy1+ii)+bc(iyy2+ii)                                  2d29s20
          bc(ixx1+ii)=fo*(xsum+ysum)
          bc(ixx2+ii)=fo*(xsum-ysum)                                         2d29s20
         end do
         ibcoff=iyy1                                                    2d29s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixz1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixz2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(ixz1+ii)=fo*(bc(ixz1+ii)+bc(ixz2+ii))                               2d26s20
         end do                                                          2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyz1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=iyz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyz2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(iyz1+ii)=fo*(bc(iyz1+ii)+bc(iyz2+ii))                               2d26s20
         end do                                                          2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixy1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixy1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        ibc(ja),bc(ja2),bc(ja3),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        lb8,bc(jb2),bc(isnormb),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixy2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(ixy1+ii)=fo*(bc(ixy1+ii)+bc(ixy2+ii))                               2d26s20
         end do                                                          2d26s20
         jzz1=izz1                                                       2d26s20
         jxx1=ixx1                                                       2d26s20
         jxx2=ixx2                                                       2d26s20
         jxz1=ixz1                                                       2d26s20
         jyz1=iyz1                                                       2d26s20
         jxy1=ixy1                                                       2d26s20
         do id=0,ndl-1                                                    2d24s20
          idp=id+iblstor(ipack2(2))                                      2d26s20
          idpp=idp+nbaslarge
          do ic=0,ncs-1                                                   2d24s20
           icp=ic+ibsstor(ipack2(1))+nbaslarge2                         2d29s20
           icpp=icp+nbassmall                                            2d26s20
           do ib=0,nbs-1                                                  2d24s20
            ibp=ib+2*nbl                                                2d27s20
            ibpp=ibp+nbs                                                2d27s20
            do ia=0,nal-1                                                 2d24s20
             iap=ia
             iapp=iap+nal                                                 2d27s20
             iad11=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))              2d27s20
             iad22=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))              2d27s20
             iad33=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))              2d27s20
             iad44=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))              2d27s20
             bc(iad11)=bc(iad11)+bc(jzz1)
             bc(iad22)=bc(iad22)-bc(jzz1)
             bc(iad33)=bc(iad33)-bc(jzz1)
             bc(iad44)=bc(iad44)+bc(jzz1)
             iad41=itmpe+iapp+naa*(ibp+nbb*(icpp+nbtot*idp))            2d27s20
             iad32=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idp))            2d27s20
             iad23=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idpp))            2d27s20
             iad14=itmpe+iap+naa*(ibpp+nbb*(icp+nbtot*idpp))            2d27s20
             bc(iad41)=bc(jxx2)                                         3d1s20
             bc(iad32)=bc(iad32)+bc(jxx1)                               3d1s20
             bc(iad23)=bc(iad23)+bc(jxx1)                               3d1s20
             bc(iad14)=bc(jxx2)                                         3d1s20
             iad41=iad41+nwds                                           2d27s20
             iad14=iad14+nwds                                           2d27s20
             bc(iad41)=2d0*bc(jxy1)                                     2d27s20
             bc(iad14)=-2d0*bc(jxy1)                                    2d27s20
             iad21=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idp))             2d27s20
             iad31=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idp))             2d27s20
             iad12=itmpe+iap+naa*(ibpp+nbb*(icp+nbtot*idp))             2d27s20
             iad13=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idpp))             2d27s20
             bc(iad21)=bc(jxz1)                                         2d27s20
             bc(iad31)=bc(jxz1)                                         2d27s20
             bc(iad12)=bc(jxz1)                                         2d27s20
             bc(iad13)=bc(jxz1)                                         2d27s20
             iad21=iad21+nwds                                           2d27s20
             iad31=iad31+nwds                                           2d27s20
             iad12=iad12+nwds                                           2d27s20
             iad13=iad13+nwds                                           2d27s20
             bc(iad21)=bc(jyz1)                                         2d27s20
             bc(iad31)=bc(jyz1)                                         2d27s20
             bc(iad12)=-bc(jyz1)                                         2d27s20
             bc(iad13)=-bc(jyz1)                                         2d27s20
             iad42=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idp))           2d27s20
             iad43=itmpe+iapp+naa*(ibp+nbb*(icpp+nbtot*idpp))           2d27s20
             iad24=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idpp))           2d27s20
             iad34=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idpp))           2d27s20
             bc(iad42)=-bc(jxz1)                                        2d27s20
             bc(iad43)=-bc(jxz1)                                        2d27s20
             bc(iad24)=-bc(jxz1)                                        2d27s20
             bc(iad34)=-bc(jxz1)                                        2d27s20
             iad42=iad42+nwds                                           2d27s20
             iad43=iad43+nwds                                           2d27s20
             iad24=iad24+nwds                                           2d27s20
             iad34=iad34+nwds                                           2d27s20
             bc(iad42)=-bc(jyz1)                                        2d27s20
             bc(iad43)=-bc(jyz1)                                        2d27s20
             bc(iad24)=bc(jyz1)                                         2d27s20
             bc(iad34)=bc(jyz1)                                         2d27s20
             jzz1=jzz1+1                                                    2d24s20
             jxx1=jxx1+1                                                    2d24s20
             jxx2=jxx2+1                                                3d1s20
             jxz1=jxz1+1                                                    2d24s20
             jyz1=jyz1+1                                                    2d24s20
             jxy1=jxy1+1                                                    2d24s20
            end do
           end do
          end do
         end do
         ibcoff=izz1                                                    3d1s20
         nla=2*la8+1                                                     2d25s20
         nlb=2*ibc(jb)+1                                                 2d26s20
         nlc=2*ibc(jc)+1                                                  2d24s20
         nld=2*ld8+1                                                     2d25s20
         xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $        *dfloat(irtyp(5))                                          2d26s20
         xnorm=-xnorm*0.5d0*scaleb                                              2d26s20
         nabcd=nas*nbl*ncl*nds                                          3d1s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,izz1,0,0,1,0,0,0,0,0,0,0,0,0,                    2d26s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=izz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,izz2,0,0,0,0,0,1,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(izz1+ii)=fo*(bc(izz1+ii)+bc(izz2+ii))
         end do                                                          2d26s20
         ibcoff=izz2                                                     2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixx1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        2,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixx1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixx2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        2,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixx2+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyy1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=iyy1+nabcd
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyy2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          xsum=bc(ixx1+ii)+bc(ixx2+ii)                                  2d29s20
          ysum=bc(iyy1+ii)+bc(iyy2+ii)                                  2d29s20
          bc(ixx1+ii)=fo*(xsum+ysum)
          bc(ixx2+ii)=fo*(xsum-ysum)                                         2d29s20
         end do
         ibcoff=iyy1                                                    2d29s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixz1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixz2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(ixz1+ii)=fo*(bc(ixz1+ii)+bc(ixz2+ii))                               2d26s20
         end do                                                          2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyz1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=iyz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyz2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(iyz1+ii)=fo*(bc(iyz1+ii)+bc(iyz2+ii))                               2d26s20
         end do                                                          2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixy1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixy1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        ibc(jc),bc(jc2),bc(jc3),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ld8,bc(jd2),bc(isnormd),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixy2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(ixy1+ii)=fo*(bc(ixy1+ii)+bc(ixy2+ii))                               2d26s20
         end do                                                          2d26s20
         jzz1=izz1                                                       2d26s20
         jxx1=ixx1                                                       2d26s20
         jxx2=ixx2                                                      3d1s20
         jxz1=ixz1                                                       2d26s20
         jyz1=iyz1                                                       2d26s20
         jxy1=ixy1                                                       2d26s20
         do id=0,nds-1                                                    2d24s20
          idp=id+ibsstor(ipack2(2))+nbaslarge2                          2d29s20
          idpp=idp+nbassmall                                             2d25s20
          do ic=0,ncl-1                                                   2d24s20
           icp=ic+iblstor(ipack2(1))                                      2d24s20
           icpp=icp+nbaslarge                                             2d24s20
           do ib=0,nbl-1                                                  2d24s20
            ibp=ib
            ibpp=ibp+nbl                                                  2d27s20
            do ia=0,nas-1                                                 2d24s20
             iap=ia+nal*2                                               2d27s20
             iapp=iap+nas                                               2d27s20
             iad11=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))              2d27s20
             iad22=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))              2d27s20
             iad33=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))              2d27s20
             iad44=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))              2d27s20
             bc(iad11)=bc(iad11)+bc(jzz1)                               3d1s20
             bc(iad22)=bc(iad22)-bc(jzz1)                               3d1s20
             bc(iad33)=bc(iad33)-bc(jzz1)                               3d1s20
             bc(iad44)=bc(iad44)+bc(jzz1)                               3d1s20
             iad41=itmpe+iapp+naa*(ibp+nbb*(icpp+nbtot*idp))            2d27s20
             iad32=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idp))            2d27s20
             iad23=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idpp))            2d27s20
             iad14=itmpe+iap+naa*(ibpp+nbb*(icp+nbtot*idpp))            2d27s20
             bc(iad41)=bc(jxx2)                                         3d1s20
             bc(iad32)=bc(iad32)+bc(jxx1)                               3d1s20
             bc(iad23)=bc(iad23)+bc(jxx1)                               3d1s20
             bc(iad14)=bc(jxx2)                                         3d1s20
             iad41=iad41+nwds                                           2d27s20
             iad14=iad14+nwds                                           2d27s20
             bc(iad41)=2d0*bc(jxy1)                                     2d27s20
             bc(iad14)=-2d0*bc(jxy1)                                    2d27s20
             iad21=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idp))             2d27s20
             iad31=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idp))             2d27s20
             iad12=itmpe+iap+naa*(ibpp+nbb*(icp+nbtot*idp))             2d27s20
             iad13=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idpp))             2d27s20
             bc(iad21)=bc(jxz1)                                         2d27s20
             bc(iad31)=bc(jxz1)                                         2d27s20
             bc(iad12)=bc(jxz1)                                         2d27s20
             bc(iad13)=bc(jxz1)                                         2d27s20
             iad21=iad21+nwds                                           2d27s20
             iad31=iad31+nwds                                           2d27s20
             iad12=iad12+nwds                                           2d27s20
             iad13=iad13+nwds                                           2d27s20
             bc(iad21)=bc(jyz1)                                         2d27s20
             bc(iad31)=bc(jyz1)                                         2d27s20
             bc(iad12)=-bc(jyz1)                                         2d27s20
             bc(iad13)=-bc(jyz1)                                         2d27s20
             iad42=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idp))           2d27s20
             iad43=itmpe+iapp+naa*(ibp+nbb*(icpp+nbtot*idpp))           2d27s20
             iad24=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idpp))           2d27s20
             iad34=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idpp))           2d27s20
             bc(iad42)=-bc(jxz1)                                        2d27s20
             bc(iad43)=-bc(jxz1)                                        2d27s20
             bc(iad24)=-bc(jxz1)                                        2d27s20
             bc(iad34)=-bc(jxz1)                                        2d27s20
             iad42=iad42+nwds                                           2d27s20
             iad43=iad43+nwds                                           2d27s20
             iad24=iad24+nwds                                           2d27s20
             iad34=iad34+nwds                                           2d27s20
             bc(iad42)=-bc(jyz1)                                        2d27s20
             bc(iad43)=-bc(jyz1)                                        2d27s20
             bc(iad24)=bc(jyz1)                                         2d27s20
             bc(iad34)=bc(jyz1)                                         2d27s20
             jzz1=jzz1+1                                                    2d24s20
             jxx1=jxx1+1                                                    2d24s20
             jxx2=jxx2+1                                                3d1s20
             jxz1=jxz1+1                                                    2d24s20
             jyz1=jyz1+1                                                    2d24s20
             jxy1=jxy1+1                                                    2d24s20
            end do
           end do
          end do
         end do
         ibcoff=izz1                                                    3d1s20
         nla=2*la8+1                                                     2d25s20
         nlb=2*ibc(jb)+1                                                 2d26s20
         nlc=2*lc8+1                                                     2d26s20
         nld=2*ibc(jd)+1                                                 2d26s20
         xnorm=oneover4pi*oneover4pi*sqrt(dfloat(nla*nlb*nlc*nld))        2d24s20
     $        *dfloat(irtyp(5))                                          2d26s20
         xnorm=-xnorm*0.5d0*scaleb                                              2d26s20
         nabcd=nas*nbl*ncs*ndl                                          3d1s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,izz1,0,0,1,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=izz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,izz2,0,0,0,0,0,1,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(izz1+ii)=fo*(bc(izz1+ii)+bc(izz2+ii))
         end do                                                          2d26s20
         ibcoff=izz2                                                     2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixx1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        2,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixx1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixx2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        2,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixx2+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyy1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=iyy1+nabcd
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyy2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          xsum=bc(ixx1+ii)+bc(ixx2+ii)                                  2d29s20
          ysum=bc(iyy1+ii)+bc(iyy2+ii)                                  2d29s20
          bc(ixx1+ii)=fo*(xsum+ysum)
          bc(ixx2+ii)=fo*(xsum-ysum)                                         2d29s20
         end do
         ibcoff=iyy1                                                    2d29s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixz1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixz2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(ixz1+ii)=fo*(bc(ixz1+ii)+bc(ixz2+ii))                               2d26s20
         end do                                                          2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyz1,0,1,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=iyz1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,iyz2,0,0,0,0,1,0,0,0,0,0,0,0,                     2d24s20
     $        4,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(iyz1+ii)=fo*(bc(iyz1+ii)+bc(iyz2+ii))                               2d26s20
         end do                                                          2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixy1,1,0,0,0,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         ibcoff=ixy1+nabcd                                               2d26s20
         call derid(                                                      2d24s20
     $        la8,bc(ja2),bc(isnorma),bc(ja5),bc(ja6),bc(ja7),ibc(ja4),   2d25s20
     $        ibc(jb),bc(jb2),bc(jb3),bc(jb5),bc(jb6),bc(jb7),ibc(jb4),   2d24s20
     $        lc8,bc(jc2),bc(isnormc),bc(jc5),bc(jc6),bc(jc7),ibc(jc4),   2d24s20
     $        ibc(jd),bc(jd2),bc(jd3),bc(jd5),bc(jd6),bc(jd7),ibc(jd4),   2d24s20
     $        dum,1,idum,ixy2,0,0,0,1,0,0,0,0,0,0,0,0,                     2d24s20
     $        3,.false.,xnorm,bc,ibc)                                   11d17s22
         do ii=0,nabcd-1                                                  2d26s20
          bc(ixy1+ii)=fo*(bc(ixy1+ii)+bc(ixy2+ii))                               2d26s20
         end do                                                          2d26s20
         jzz1=izz1                                                       2d26s20
         jxx1=ixx1                                                       2d26s20
         jxx2=ixx2                                                      3d1s20
         jxz1=ixz1                                                       2d26s20
         jyz1=iyz1                                                       2d26s20
         jxy1=ixy1                                                       2d26s20
         do id=0,ndl-1                                                    2d24s20
          idp=id+iblstor(ipack2(2))                                      2d26s20
          idpp=idp+nbaslarge                                             2d26s20
          do ic=0,ncs-1                                                   2d24s20
           icp=ic+ibsstor(ipack2(1))+nbaslarge2                         2d29s20
           icpp=icp+nbassmall                                            2d26s20
           do ib=0,nbl-1                                                  2d24s20
            ibp=ib
            ibpp=ibp+nbl                                                  2d27s20
            do ia=0,nas-1                                                 2d24s20
             iap=ia+nal*2                                               2d27s20
             iapp=iap+nas                                               2d27s20
             iad11=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idp))              2d27s20
             iad22=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idp))              2d27s20
             iad33=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idpp))              2d27s20
             iad44=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idpp))              2d27s20
             bc(iad11)=bc(iad11)+bc(jzz1)                               3d1s20
             bc(iad22)=bc(iad22)-bc(jzz1)                               3d1s20
             bc(iad33)=bc(iad33)-bc(jzz1)                               3d1s20
             bc(iad44)=bc(iad44)+bc(jzz1)                               3d1s20
             iad41=itmpe+iapp+naa*(ibp+nbb*(icpp+nbtot*idp))            2d27s20
             iad32=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idp))            2d27s20
             iad23=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idpp))            2d27s20
             iad14=itmpe+iap+naa*(ibpp+nbb*(icp+nbtot*idpp))            2d27s20
             bc(iad41)=bc(jxx2)                                         3d1s20
             bc(iad32)=bc(iad32)+bc(jxx1)                               3d1s20
             bc(iad23)=bc(iad23)+bc(jxx1)                               3d1s20
             bc(iad14)=bc(jxx2)                                         3d1s20
             iad41=iad41+nwds                                           2d27s20
             iad14=iad14+nwds                                           2d27s20
             bc(iad41)=2d0*bc(jxy1)                                     2d27s20
             bc(iad14)=-2d0*bc(jxy1)                                    2d27s20
             iad21=itmpe+iapp+naa*(ibp+nbb*(icp+nbtot*idp))             2d27s20
             iad31=itmpe+iap+naa*(ibp+nbb*(icpp+nbtot*idp))             2d27s20
             iad12=itmpe+iap+naa*(ibpp+nbb*(icp+nbtot*idp))             2d27s20
             iad13=itmpe+iap+naa*(ibp+nbb*(icp+nbtot*idpp))             2d27s20
             bc(iad21)=bc(jxz1)                                         2d27s20
             bc(iad31)=bc(jxz1)                                         2d27s20
             bc(iad12)=bc(jxz1)                                         2d27s20
             bc(iad13)=bc(jxz1)                                         2d27s20
             iad21=iad21+nwds                                           2d27s20
             iad31=iad31+nwds                                           2d27s20
             iad12=iad12+nwds                                           2d27s20
             iad13=iad13+nwds                                           2d27s20
             bc(iad21)=bc(jyz1)                                         2d27s20
             bc(iad31)=bc(jyz1)                                         2d27s20
             bc(iad12)=-bc(jyz1)                                         2d27s20
             bc(iad13)=-bc(jyz1)                                         2d27s20
             iad42=itmpe+iapp+naa*(ibpp+nbb*(icpp+nbtot*idp))           2d27s20
             iad43=itmpe+iapp+naa*(ibp+nbb*(icpp+nbtot*idpp))           2d27s20
             iad24=itmpe+iapp+naa*(ibpp+nbb*(icp+nbtot*idpp))           2d27s20
             iad34=itmpe+iap+naa*(ibpp+nbb*(icpp+nbtot*idpp))           2d27s20
             bc(iad42)=-bc(jxz1)                                        2d27s20
             bc(iad43)=-bc(jxz1)                                        2d27s20
             bc(iad24)=-bc(jxz1)                                        2d27s20
             bc(iad34)=-bc(jxz1)                                        2d27s20
             iad42=iad42+nwds                                           2d27s20
             iad43=iad43+nwds                                           2d27s20
             iad24=iad24+nwds                                           2d27s20
             iad34=iad34+nwds                                           2d27s20
             bc(iad42)=-bc(jyz1)                                        2d27s20
             bc(iad43)=-bc(jyz1)                                        2d27s20
             bc(iad24)=bc(jyz1)                                         2d27s20
             bc(iad34)=bc(jyz1)                                         2d27s20
             jzz1=jzz1+1                                                    2d24s20
             jxx1=jxx1+1                                                    2d24s20
             jxx2=jxx2+1                                                3d1s20
             jxz1=jxz1+1                                                    2d24s20
             jyz1=jyz1+1                                                    2d24s20
             jxy1=jxy1+1                                                    2d24s20
            end do
           end do
          end do
         end do
         ibcoff=izz1
        end if                                                           2d26s20
       end do                                                           2d27s20
c
c     we have all integrals for a,b so
c     transform c&d
c
       naabb=naa*nbb                                                     2d28s20
       nrow=naabb*nbtot                                                  2d28s20
       nwds=nrow*nox                                                     2d28s20
       itmpr=ibcoff                                                      2d28s20
       itmpi=itmpr+nwds                                                  2d28s20
       ibcoff=itmpi+nwds                                                 2d28s20
       call enough('paraeri4ct.  4',bc,ibc)
       call dgemm('n','n',nrow,nox,nbtot,1d0,bc(itmpe),nrow,vec,nbtot,   2d28s20
     $      0d0,bc(itmpr),nrow,                                          2d28s20
     d' paraeri4ct.  1')
       call dgemm('n','n',nrow,nox,nbtot,1d0,bc(itmpe),nrow,vec(1,1,2),  2d28s20
     $      nbtot,0d0,bc(itmpi),nrow,                                    2d28s20
     d' paraeri4ct.  2')
       if(irtyp(5).ne.0)then                                             2d28s20
        call dgemm('n','n',nrow,nox,nbtot,-1d0,bc(itmpei),nrow,          2d28s20
     $      vec(1,1,2),nbtot,1d0,bc(itmpr),nrow,                         2d28s20
     d' paraeri4ct.  3')
        call dgemm('n','n',nrow,nox,nbtot,1d0,bc(itmpei),nrow,vec,       2d28s20
     $       nbtot,1d0,bc(itmpi),nrow,                                    2d28s20
     d' paraeri4ct.  4')
       end if                                                            2d28s20
       itmp2r=ibcoff                                                     2d28s20
       itmp2i=itmp2r+nwds                                                2d28s20
       ibcoff=itmp2i+nwds                                                2d28s20
       do ii=0,nox-1                                                      2d28s20
        do j=0,nbtot-1                                                   2d28s20
         jir=itmpr+naabb*(j+nbtot*ii)                                      2d28s20
         jii=jir+nwds                                                    2d28s20
         ijr=itmp2r+naabb*(ii+nox*j)                                      2d28s20
         iji=ijr+nwds                                                    2d28s20
         do k=0,naabb-1                                                  2d28s20
          bc(ijr+k)=bc(jir+k)                                            2d28s20
          bc(iji+k)=+bc(jii+k)                                          2d28s20
         end do                                                          2d28s20
        end do                                                           2d28s20
       end do                                                            2d28s20
       nrow=naabb*nox                                                    2d28s20
       nwds=nrow*nox                                                    2d28s20
       itmpi=itmpr+nwds                                                  2d28s20
       call dgemm('n','n',nrow,nox,nbtot,1d0,bc(itmp2r),nrow,vec,nbtot,  2d28s20
     $      0d0,bc(itmpr),nrow,                                          2d28s20
     d' paraeri4ct.  5')
       call dgemm('n','n',nrow,nox,nbtot,+1d0,bc(itmp2i),nrow,          2d28s20
     $      vec(1,1,2),nbtot,1d0,bc(itmpr),nrow,                        2d28s20
     d' paraeri4ct.  6')
       call dgemm('n','n',nrow,nox,nbtot,1d0,bc(itmp2i),nrow,vec,nbtot,  2d28s20
     $      0d0,bc(itmpi),nrow,                                          2d28s20
     d' paraeri4ct.  7')
       call dgemm('n','n',nrow,nox,nbtot,-1d0,bc(itmp2r),nrow,          2d28s20
     $      vec(1,1,2),nbtot,1d0,bc(itmpi),nrow,                        2d28s20
     d' paraeri4ct.  8')
       do id=0,nox-1                                                    2d28s20
        do ic=0,nox-1                                                   2d28s20
         jtmpr=itmpr+naabb*(id+nox*ic)                                  2d28s20
         do ib=0,nbl-1                                                  2d28s20
          ibp=ib+iblstor(ishellb)                                       2d28s20
          ibpp=ibp+nbaslarge                                            2d28s20
          lb=ib+nbl                                                     2d28s20
          do ia=0,nal-1                                                 2d28s20
           iap=ia+iblstor(ishella)                                      2d28s20
           iapp=iap+nbaslarge                                           2d28s20
           la=ia+nal                                                    2d28s20
           k=jtmpr+ia+naa*ib                                            2d28s20
           iad1=ihalfr+ic+nox*(id+nox*(iap+nbtot*ibp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+ia+naa*lb
           iad1=ihalfr+ic+nox*(id+nox*(iap+nbtot*ibpp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+la+naa*ib
           iad1=ihalfr+ic+nox*(id+nox*(iapp+nbtot*ibp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+la+naa*lb
           iad1=ihalfr+ic+nox*(id+nox*(iapp+nbtot*ibpp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
          end do                                                        2d28s20
          do ia=0,nas-1                                                 2d28s20
           iap=ia+ibsstor(ishella)+nbaslarge2                           2d29s20
           iapp=iap+nbassmall                                           2d28s20
           lap=ia+nal*2                                                 2d28s20
           lapp=lap+nas                                                 2d28s20
           k=jtmpr+lap+naa*ib                                            2d28s20
           iad1=ihalfr+ic+nox*(id+nox*(iap+nbtot*ibp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+lap+naa*lb
           iad1=ihalfr+ic+nox*(id+nox*(iap+nbtot*ibpp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+lapp+naa*ib
           iad1=ihalfr+ic+nox*(id+nox*(iapp+nbtot*ibp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+lapp+naa*lb
           iad1=ihalfr+ic+nox*(id+nox*(iapp+nbtot*ibpp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
          end do                                                        2d28s20
         end do                                                         2d28s20
         do ib=0,nbs-1                                                  2d28s20
          ibp=ib+ibsstor(ishellb)+nbaslarge2                            2d29s20
          ibpp=ibp+nbassmall                                            2d28s20
          lbp=ib+2*nbl                                                  2d28s20
          lbpp=lbp+nbs                                                  2d28s20
          do ia=0,nal-1                                                 2d28s20
           iap=ia+iblstor(ishella)                                      2d28s20
           iapp=iap+nbaslarge                                           2d28s20
           la=ia+nal                                                    2d28s20
           k=jtmpr+ia+naa*lbp                                            2d28s20
           iad1=ihalfr+ic+nox*(id+nox*(iap+nbtot*ibp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+ia+naa*lbpp
           iad1=ihalfr+ic+nox*(id+nox*(iap+nbtot*ibpp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+la+naa*lbp
           iad1=ihalfr+ic+nox*(id+nox*(iapp+nbtot*ibp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+la+naa*lbpp
           iad1=ihalfr+ic+nox*(id+nox*(iapp+nbtot*ibpp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
          end do                                                        2d28s20
          do ia=0,nas-1                                                 2d28s20
           iap=ia+ibsstor(ishella)+nbaslarge2                           2d29s20
           iapp=iap+nbassmall                                           2d28s20
           lap=ia+nal*2                                                 2d28s20
           lapp=lap+nas                                                 2d28s20
           k=jtmpr+lap+naa*lbp                                            2d28s20
           iad1=ihalfr+ic+nox*(id+nox*(iap+nbtot*ibp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+lap+naa*lbpp
           iad1=ihalfr+ic+nox*(id+nox*(iap+nbtot*ibpp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+lapp+naa*lbp
           iad1=ihalfr+ic+nox*(id+nox*(iapp+nbtot*ibp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
           k=jtmpr+lapp+naa*lbpp
           iad1=ihalfr+ic+nox*(id+nox*(iapp+nbtot*ibpp))                  2d28s20
           bc(iad1)=bc(k)                                               2d28s20
           k=k+nwds                                                     2d28s20
           iad1=iad1+nwdsh
           bc(iad1)=+bc(k)                                               2d28s20
          end do                                                        2d28s20
         end do                                                         2d28s20
        end do                                                          2d28s20
       end do                                                           2d28s20
       ibcoff=itmpe                                                     2d28s20
      end do                                                            2d19s10
      call dws_gsumf(bc(ihalfr),nwdsh*2)                                2d28s20
      rmsab=0d0
      rmsabi=0d0
      nok=0
      do ia=0,nbtot-1
       do ib=0,ia-1
        do ic=0,nox-1
         do id=0,ic-1
          iad11=ic+nox*(id+nox*(ia+nbtot*ib))
          iad21=id+nox*(ic+nox*(ib+nbtot*ia))
          rmsab=rmsab+(bc(ihalfr+iad11)-bc(ihalfr+iad21))**2
          rmsabi=rmsabi+(bc(ihalfi+iad11)+bc(ihalfi+iad21))**2
         end do
        end do
       end do
      end do
      rmsab=sqrt(rmsab/(nbtot*nbtot*nox*nox))
      rmsabi=sqrt(rmsabi/(nbtot*nbtot*nox*nox))
c     now transform b
c
      itmpr=ibcoff                                                      2d28s20
      nxx=nox*nox                                                       2d28s20
      nrow=nxx*nbtot                                                    2d28s20
      nwds=nrow*nox                                                     2d28s20
      itmpi=itmpr+nwds                                                  2d28s20
      ibcoff=itmpi+nwds                                                 2d28s20
      call enough('paraeri4ct.  5',bc,ibc)
      call dgemm('n','n',nrow,nox,nbtot,1d0,bc(ihalfr),nrow,vec,nbtot,  2d28s20
     $     0d0,bc(itmpr),nrow,                                          2d28s20
     d' paraeri4ct.  9')
      call dgemm('n','n',nrow,nox,nbtot,-1d0,bc(ihalfi),nrow,           2d28s20
     $     vec(1,1,2),nbtot,1d0,bc(itmpr),nrow,                         2d28s20
     d' paraeri4ct. 10')
      call dgemm('n','n',nrow,nox,nbtot,1d0,bc(ihalfi),nrow,vec,nbtot,  2d28s20
     $     0d0,bc(itmpi),nrow,                                          2d28s20
     d' paraeri4ct. 11')
      call dgemm('n','n',nrow,nox,nbtot,1d0,bc(ihalfr),nrow,vec(1,1,2), 2d28s20
     $     nbtot,1d0,bc(itmpi),nrow,                                    2d28s20
     d' paraeri4ct. 12')
      do i=0,nox-1                                                      2d28s20
       do j=0,nbtot-1                                                   2d28s20
        jir=itmpr+nxx*(j+nbtot*i)                                       2d28s20
        jii=jir+nwds                                                    2d28s20
        ijr=ihalfr+nxx*(i+nox*j)                                        2d28s20
        iji=ijr+nwds                                                    2d28s20
        do k=0,nxx-1                                                    2d28s20
         bc(ijr+k)=bc(jir+k)                                            2d28s20
         bc(iji+k)=+bc(jii+k)                                           2d28s20
        end do                                                          2d28s20
       end do                                                           2d28s20
      end do                                                            2d28s20
      ihalfi=ihalfr+nwds                                                2d28s20
      nxx2=nxx*nxx                                                      2d28s20
      itmpi=itmpr+nxx2                                                  2d28s20
      nrow=nxx*nox                                                      2d28s20
      call dgemm('n','n',nrow,nox,nbtot,1d0,bc(ihalfr),nrow,vec,nbtot,  2d28s20
     $     0d0,bc(itmpr),nrow,                                          2d28s20
     d' paraeri4ct. 13')
      call dgemm('n','n',nrow,nox,nbtot,+1d0,bc(ihalfi),nrow,vec(1,1,2),2d28s20
     $     nbtot,1d0,bc(itmpr),nrow,                                    2d28s20
     d' paraeri4ct. 14')
      call dgemm('n','n',nrow,nox,nbtot,-1d0,bc(ihalfr),nrow,vec(1,1,2),2d28s20
     $     nbtot,0d0,bc(itmpi),nrow,                                    2d28s20
     d' paraeri4ct. 15')
      call dgemm('n','n',nrow,nox,nbtot,1d0,bc(ihalfi),nrow,vec,        2d28s20
     $     nbtot,1d0,bc(itmpi),nrow,                                    2d28s20
     d' paraeri4ct. 16')
      do ia=0,nox-1                                                      2d28s20
       do ib=0,nox-1                                                     2d28s20
        do k=0,nxx-1                                                    2d28s20
         iad1=itmpr+k+nxx*(ib+nox*ia)                                   2d28s20
         iad2=ia+1+nox*(ib+nox*k)                                       2d28s20
         oooo(iad2)=bc(iad1)                                            2d28s20
         iad1=iad1+nxx2
         iad2=iad2+nxx2                                                 2d28s20
         oooo(iad2)=+bc(iad1)                                           2d28s20
        end do                                                          2d28s20
       end do                                                           2d28s20
      end do                                                             2d28s20
      ibcoff=ibcoffo                                                    2d19s10
      return
      end                                                               2d19s10
