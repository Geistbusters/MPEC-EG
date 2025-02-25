c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine restid(la,zza,xna,xa,ya,za,iadda,                         1d29s10
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               icartt,iax,iay,iaz,ibx,iby,ibz,icx,icy,icz,        4d27s16
     $     idx,idy,idz,fmulq,idorel,ascale,ldeb,bc,ibc)                 11d10s22
      implicit real*8 (a-h,o-z)
      logical ldeb,lab,lcd                                              10d8s15
c
c     take electron repulsion integrals and store them properly
c     for 2 component calculation, and optionally compute additional
c     2 electron integrals.
c     this version includes derivatives for geometry dependence of      4d27s16
c     quantities.                                                       4d27s16
c
      include "common.store"
      include "common.spher"
      include "common.rys"                                              7d5s12
      dimension iovrm(8),ikinm(8),ioffb(8),ioffk(8),nhere(8,4),         1d13s10
     $     mhere(8,4),icoef(8,4),ipow(8,4),ltmp(4),mult(8,8),joff(8,4)  1d14s10
      dimension eraw(nbasall,nbasall,1),nxyz(3,3),breit(12),ibt(13,12,4)10d8s15
      common/timerocm/tovr,telapo(15)                                    5d7s12
      data icall/0/                                                     7d22s14
      data mult/1,2,3,4,5,6,7,8, 2,1,4,3,6,5,8,7,                       1d13s10
     $          3,4,1,2,7,8,5,6, 4,3,2,1,8,7,6,5,                       1d13s10
     $          5,6,7,8,1,2,3,4, 6,5,8,7,2,1,4,3,                       1d13s10
     $          7,8,5,6,3,4,1,2, 8,7,6,5,4,3,2,1/                       1d13s10
      data nxyz/1,0,0, 0,1,0, 0,0,1/                                    10d8s15
      data breit/3*1d0,9*-1d0/                                          10d8s15
      data (ibt(j,1,1),j=1,13) /1,0,0, 0,0,0, 1,0,0, 0,0,0, 1/
      data (ibt(j,2,1),j=1,13) /0,1,0, 0,0,0, 0,1,0, 0,0,0, 1/
      data (ibt(j,3,1),j=1,13) /0,0,1, 0,0,0, 0,0,1, 0,0,0, 1/
      data (ibt(j,4,1),j=1,13) /1,0,0, 0,0,0, 0,0,1, 0,0,0, 6/
      data (ibt(j,5,1),j=1,13) /1,0,0, 0,0,0, 1,0,0, 0,0,0, 2/
      data (ibt(j,6,1),j=1,13) /1,0,0, 0,0,0, 0,1,0, 0,0,0, 5/
      data (ibt(j,7,1),j=1,13) /0,1,0, 0,0,0, 1,0,0, 0,0,0, 5/
      data (ibt(j,8,1),j=1,13) /0,1,0, 0,0,0, 0,1,0, 0,0,0, 3/
      data (ibt(j,9,1),j=1,13) /0,1,0, 0,0,0, 0,0,1, 0,0,0, 7/
      data (ibt(j,10,1),j=1,13)/0,0,1, 0,0,0, 1,0,0, 0,0,0, 6/
      data (ibt(j,11,1),j=1,13)/0,0,1, 0,0,0, 0,1,0, 0,0,0, 7/
      data (ibt(j,12,1),j=1,13)/0,0,1, 0,0,0, 0,0,1, 0,0,0, 4/
      data (ibt(j,1,2),j=1,13) /0,0,0, 1,0,0, 0,0,0, 1,0,0, 1/
      data (ibt(j,2,2),j=1,13) /0,0,0, 0,1,0, 0,0,0, 0,1,0, 1/
      data (ibt(j,3,2),j=1,13) /0,0,0, 0,0,1, 0,0,0, 0,0,1, 1/
      data (ibt(j,4,2),j=1,13) /0,0,0, 1,0,0, 0,0,0, 0,0,1, 6/
      data (ibt(j,5,2),j=1,13) /0,0,0, 1,0,0, 0,0,0, 1,0,0, 2/
      data (ibt(j,6,2),j=1,13) /0,0,0, 1,0,0, 0,0,0, 0,1,0, 5/
      data (ibt(j,7,2),j=1,13) /0,0,0, 0,1,0, 0,0,0, 1,0,0, 5/
      data (ibt(j,8,2),j=1,13) /0,0,0, 0,1,0, 0,0,0, 0,1,0, 3/
      data (ibt(j,9,2),j=1,13) /0,0,0, 0,1,0, 0,0,0, 0,0,1, 7/
      data (ibt(j,10,2),j=1,13)/0,0,0, 0,0,1, 0,0,0, 1,0,0, 6/
      data (ibt(j,11,2),j=1,13)/0,0,0, 0,0,1, 0,0,0, 0,1,0, 7/
      data (ibt(j,12,2),j=1,13)/0,0,0, 0,0,1, 0,0,0, 0,0,1, 4/
      data (ibt(j,1,3),j=1,13) /0,0,0, 1,0,0, 1,0,0, 0,0,0, 1/
      data (ibt(j,2,3),j=1,13) /0,0,0, 0,1,0, 0,1,0, 0,0,0, 1/
      data (ibt(j,3,3),j=1,13) /0,0,0, 0,0,1, 0,0,1, 0,0,0, 1/
      data (ibt(j,4,3),j=1,13) /0,0,0, 1,0,0, 0,0,1, 0,0,0, 6/
      data (ibt(j,5,3),j=1,13) /0,0,0, 1,0,0, 1,0,0, 0,0,0, 2/
      data (ibt(j,6,3),j=1,13) /0,0,0, 1,0,0, 0,1,0, 0,0,0, 5/
      data (ibt(j,7,3),j=1,13) /0,0,0, 0,1,0, 1,0,0, 0,0,0, 5/
      data (ibt(j,8,3),j=1,13) /0,0,0, 0,1,0, 0,1,0, 0,0,0, 3/
      data (ibt(j,9,3),j=1,13) /0,0,0, 0,1,0, 0,0,1, 0,0,0, 7/
      data (ibt(j,10,3),j=1,13)/0,0,0, 0,0,1, 1,0,0, 0,0,0, 6/
      data (ibt(j,11,3),j=1,13)/0,0,0, 0,0,1, 0,1,0, 0,0,0, 7/
      data (ibt(j,12,3),j=1,13)/0,0,0, 0,0,1, 0,0,1, 0,0,0, 4/
      data (ibt(j,1,4),j=1,13) /1,0,0, 0,0,0, 0,0,0, 1,0,0, 1/
      data (ibt(j,2,4),j=1,13) /0,1,0, 0,0,0, 0,0,0, 0,1,0, 1/
      data (ibt(j,3,4),j=1,13) /0,0,1, 0,0,0, 0,0,0, 0,0,1, 1/
      data (ibt(j,4,4),j=1,13) /1,0,0, 0,0,0, 0,0,0, 0,0,1, 6/
      data (ibt(j,5,4),j=1,13) /1,0,0, 0,0,0, 0,0,0, 1,0,0, 2/
      data (ibt(j,6,4),j=1,13) /1,0,0, 0,0,0, 0,0,0, 0,1,0, 5/
      data (ibt(j,7,4),j=1,13) /0,1,0, 0,0,0, 0,0,0, 1,0,0, 5/
      data (ibt(j,8,4),j=1,13) /0,1,0, 0,0,0, 0,0,0, 0,1,0, 3/
      data (ibt(j,9,4),j=1,13) /0,1,0, 0,0,0, 0,0,0, 0,0,1, 7/
      data (ibt(j,10,4),j=1,13)/0,0,1, 0,0,0, 0,0,0, 1,0,0, 6/
      data (ibt(j,11,4),j=1,13)/0,0,1, 0,0,0, 0,0,0, 0,1,0, 7/
      data (ibt(j,12,4),j=1,13)/0,0,1, 0,0,0, 0,0,0, 0,0,1, 4/
      save                                                              7d22s14
      ibcoffo=ibcoff                                                    8d26s15
      nsza=2*la+1
      nszac=nsza*2
      nszb=2*lb+1
      nszbc=nszb*2
      nszc=2*lc+1
      nszcc=nszc*2
      nszd=2*ld+1
      nszdc=nszd*2
c
c     first move LLLL integrals.
c
      if(nsymb.eq.-11)then
c                                                                       8d26s15
c     for no symmetry case, we are already ok for d&c, but not          8d26s15
c     for a&b.                                                          8d26s15
c                                                                       8d26s15
       do ia=nsza-1,0,-1
        do ib=nszb-1,0,-1
         ico=ib+1+nszb*ia
         icn=ib+1+nszbc*ia
         do ic=1,nbasall
          do id=1,nbasall
           eraw(id,ic,icn)=eraw(id,ic,ico)
          end do
         end do
        end do
       end do
      else                                                              8d26s15
c
c     under icartt, we have eri stored nszd,nszc,nszb,nsza.
c     move to nszdc,nszcc,nszbc,nszac.
c
       icn=icartt+nsza*nszb*nszc*nszd
       ibcoff=icn+nszac*nszbc*nszcc*nszdc
       call enough('restid.  1',bc,ibc)
       do i=0,nszac*nszbc*nszcc*nszdc-1                                 10d6s15
        bc(icn+i)=0d0                                                   10d6s15
       end do                                                           10d6s15
       jcn=icartt
       do ia=0,nsza-1
        do ib=0,nszb-1
         do ic=0,nszc-1
          do id=0,nszd-1
           iad1=icn+id+nszdc*(ic+nszcc*(ib+nszbc*ia))
           bc(iad1)=bc(jcn)
           jcn=jcn+1
          end do
         end do
        end do
       end do
       do i=0,nszac*nszbc*nszcc*nszdc-1
        bc(icartt+i)=bc(icn+i)
       end do                                                           8d26s15
       ibcoff=icartt+nszac*nszbc*nszcc*nszdc                            8d26s15
      end if                                                            8d26s15
c
c      no other 2e integrals, so zero out small component part.
c
       if(nsymb.eq.-1)then                                               10d6s15
        do ia=0,nsza-1                                                  10d6s15
         do ib=0,nszb-1                                                 10d6s15
          icsa=ib+1+nszbc*ia+nszb                                       10d6s15
          icsb=ib+1+nszbc*(ia+nsza)                                     10d6s15
          icss=ib+1+nszbc*(ia+nsza)+nszb                                10d6s15
          do ic=1,nbasall                                               10d6s15
           do id=1,nbasall                                              10d6s15
            eraw(id,ic,icsa)=0d0                                        10d6s15
            eraw(id,ic,icsb)=0d0                                        10d6s15
            eraw(id,ic,icss)=0d0                                        10d6s15
           end do                                                       10d6s15
          end do                                                        10d6s15
         end do                                                         10d6s15
        end do                                                          10d6s15
        nbh=nbasall/2
        do ia=0,nsza-1
         do ib=0,nszb-1
          icol0=ib+1+nszbc*ia
          icol1=icol0+nszb
          icol2=icol0+nszbc*nsza
          icol3=icol2+nszb
          do ic=1,nbasall
           do id=1,nbasall
            eraw(id,ic,icol1)=0d0
            eraw(id,ic,icol2)=0d0
            eraw(id,ic,icol3)=0d0
           end do
          end do
          do ic=1,nbh
           icp=ic+nbh
           do id=1,nbh
            idp=id+nbh
            eraw(idp,ic,icol0)=0d0
            eraw(idp,ic,icol1)=0d0
            eraw(idp,ic,icol2)=0d0
            eraw(id,icp,icol0)=0d0
            eraw(id,icp,icol1)=0d0
            eraw(id,icp,icol2)=0d0
            eraw(idp,icp,icol0)=0d0
           end do
          end do
         end do
        end do
       else
        do ia=0,nsza-1
         iap=ia+nsza
         do ib=0,nszb-1
          ibp=ib+nszb
          do ic=0,nszc-1
           icp=ic+nszc
           do id=0,nszd-1
            jd=id+icartt
            jdp=jd+nszd
            iad1=jd+nszdc*(ic+nszcc*(ib+nszbc*iap))
            iad2=jd+nszdc*(ic+nszcc*(ibp+nszbc*ia))
            iad3=jd+nszdc*(icp+nszcc*(ib+nszbc*ia))
            iad4=jdp+nszdc*(ic+nszcc*(ib+nszbc*ia))
            iad5=jd+nszdc*(ic+nszcc*(ibp+nszbc*iap))
            iad6=jd+nszdc*(icp+nszcc*(ib+nszbc*iap))
            iad7=jdp+nszdc*(ic+nszcc*(ib+nszbc*iap))
            iad8=jd+nszdc*(icp+nszcc*(ibp+nszbc*ia))
            iad9=jdp+nszdc*(ic+nszcc*(ibp+nszbc*ia))
            iada=jdp+nszdc*(icp+nszcc*(ib+nszbc*ia))
            iadb=jdp+nszdc*(icp+nszcc*(ibp+nszbc*ia))
            iadc=jdp+nszdc*(icp+nszcc*(ib+nszbc*iap))
            iadd=jdp+nszdc*(ic+nszcc*(ibp+nszbc*iap))
            iade=jd+nszdc*(icp+nszcc*(ibp+nszbc*iap))
            iadf=jdp+nszdc*(icp+nszcc*(ibp+nszbc*iap))
            bc(iad1)=0d0
            bc(iad2)=0d0
            bc(iad3)=0d0
            bc(iad4)=0d0
            bc(iad5)=0d0
            bc(iad6)=0d0
            bc(iad7)=0d0
            bc(iad8)=0d0
            bc(iad9)=0d0
            bc(iada)=0d0
            bc(iadb)=0d0
            bc(iadc)=0d0
            bc(iadd)=0d0
            bc(iade)=0d0
            bc(iadf)=0d0
           end do
          end do
         end do
        end do
       end if
       if(iabs(idorel).ge.2)then
c
c     ss integrals
c
       ntemp=nsza*nszb*nszc*nszd                                        10d8s15
       itemp=ibcoff                                                     10d8s15
       itemp2=itemp+ntemp                                               10d8s15
       ibcoff=itemp2+ntemp                                              10d8s15
       call enough('restid.  2',bc,ibc)
       fmul=ascale                                                      10d8s15
       icxp=icx+1                                                       4d27s16
       idxp=idx+1                                                       4d27s16
       call derid(la,zza,xna,xa,ya,za,iadda,                            10d8s15
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               itemp2,iax,iay,iax, ibx,iby,ibz, icxp,icy,icz,     4d27s16
     $               idxp,idy,idz, 1,ldeb,fmul,bc,ibc)                  11d14s22
       do i=0,ntemp-1                                                   10d8s15
        bc(itemp+i)=bc(itemp2+i)                                        10d8s15
       end do                                                           10d8s15
       icyp=icy+1                                                       4d27s16
       idyp=idy+1                                                       4d27s16
       call derid(la,zza,xna,xa,ya,za,iadda,                            10d8s15
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               itemp2,iax,iay,iaz, ibx,iby,ibz icx,icyp,icz,      4d27s16
     $               idx,idyp,idz, 1,ldeb,fmul,bc,ibc)                  11d14s22
       do i=0,ntemp-1                                                   10d8s15
        bc(itemp+i)=bc(itemp+i)+bc(itemp2+i)                            10d8s15
       end do                                                           10d8s15
       iczp=icz+1                                                       4d27s16
       idzp=idz+1                                                       4d27s16
       call derid(la,zza,xna,xa,ya,za,iadda,                            10d8s15
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               itemp2,iax,iay,iaz, ibx,iby,ibz, icx,icy,iczp,     4d27s16
     $               idx,idy,idzp, 1,ldeb,fmul,bc,ibc)                  11d14s22
       do i=0,ntemp-1                                                   10d8s15
        bc(itemp+i)=bc(itemp+i)+bc(itemp2+i)                            10d8s15
       end do                                                           10d8s15
       jtemp=itemp                                                      10d8s15
       do ia=0,nsza-1                                                   10d8s15
        do ib=0,nszb-1                                                  10d8s15
         do ic=0,nszc-1                                                 10d8s15
          icp=ic+nszc                                                   10d8s15
          do id=0,nszd-1                                                10d8s15
           jd=id+icartt                                                 10d8s15
           jdp=jd+nszd                                                  10d8s15
           iada=jdp+nszdc*(icp+nszcc*(ib+nszbc*ia))                     10d8s15
           bc(iada)=bc(jtemp)                                           10d8s15
           jtemp=jtemp+1                                                10d8s15
          end do                                                        10d8s15
         end do                                                         10d8s15
        end do                                                          10d8s15
       end do                                                           10d8s15
        iaxp=iax+1                                                      4d27s16
        ibxp=ibx+1                                                      4d27s16
        call derid(la,zza,xna,xa,ya,za,iadda,                            10d8s15
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               itemp2,iaxp,iay,iaz, ibxp,iby,ibz, icx,icy,icz,    4d27s16
     $               idx,idy,idz, 1,ldeb,fmul,bc,ibc)                   11d14s22
        do i=0,ntemp-1                                                   10d8s15
         bc(itemp+i)=bc(itemp2+i)                                        10d8s15
        end do                                                           10d8s15
        iayp=iay+1                                                      4d27s16
        ibyp=iby+1                                                      4d27s16
        call derid(la,zza,xna,xa,ya,za,iadda,                            10d8s15
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               itemp2,iax,iayp,iaz, ibx,ibyp,ibz, icx,icy,icz,     4d27s16
     $               idx,idy,idz, 1,ldeb,fmul,bc,ibc)                   11d14s22
        do i=0,ntemp-1                                                   10d8s15
         bc(itemp+i)=bc(itemp+i)+bc(itemp2+i)                            10d8s15
        end do                                                           10d8s15
        iazp=iaz+1                                                      4d27s16
        ibzp=ibz+1                                                      4d27s16
        call derid(la,zza,xna,xa,ya,za,iadda,                            10d8s15
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               itemp2,iax,iay,iazp, ibx,iby,ibzp, icx,icy,icz,    4d27s16
     $               idx,idy,idz, 1,ldeb,fmul,bc,ibc)                   11d14s22
        do i=0,ntemp-1                                                   10d8s15
         bc(itemp+i)=bc(itemp+i)+bc(itemp2+i)                            10d8s15
        end do                                                           10d8s15
        jtemp=itemp                                                      10d8s15
        do ia=0,nsza-1                                                  10d8s15
         iap=ia+nsza                                                    10d8s15
         do ib=0,nszb-1                                                 10d8s15
          ibp=ib+nszb                                                   10d8s15
          do ic=0,nszc-1                                                10d8s15
           do id=0,nszd-1                                               10d8s15
            jd=id+icartt                                                10d8s15
            iada=jd+nszdc*(ic+nszcc*(ibp+nszbc*iap))                    10d8s15
            bc(iada)=bc(jtemp)                                          10d8s15
            jtemp=jtemp+1                                               10d8s15
           end do                                                       10d8s15
          end do                                                        10d8s15
         end do                                                         10d8s15
        end do                                                          10d8s15
       ibcoff=itemp                                                     10d8s15
      end if
      if(iabs(idorel).eq.3)then
c
c     ssss integrals
c
       ntemp=nsza*nszb*nszc*nszd                                        10d8s15
       itemp=ibcoff                                                     10d8s15
       itemp2=itemp+ntemp                                               10d8s15
       ibcoff=itemp2+ntemp                                              10d8s15
       call enough('restid.  3',bc,ibc)
       fmul=ascale*ascale                                               10d8s15
       do i=0,ntemp-1                                                   10d8s15
        bc(itemp+i)=0d0                                                 10d8s15
       end do                                                           10d8s15
       do ixyz=1,3                                                      10d8s15
        do jxyz=1,3                                                     10d8s15
         iaxp=iax+nxyz(1,ixyz)
         iayp=iay+nxyz(2,ixyz)
         iazp=iaz+nxyz(3,ixyz)
         ibxp=ibx+nxyz(1,ixyz)
         ibyp=iby+nxyz(2,ixyz)
         ibzp=ibz+nxyz(3,ixyz)
         icxp=icx+nxyz(1,jxyz)
         icyp=icy+nxyz(2,jxyz)
         iczp=icz+nxyz(3,jxyz)
         idxp=idx+nxyz(1,jxyz)
         idyp=idy+nxyz(2,jxyz)
         idzp=idz+nxyz(3,jxyz)
         call derid(la,zza,xna,xa,ya,za,iadda,                            10d8s15
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               itemp2,iaxp,iayp,iazp, ibxp,ibyp,ibzp,             4d27s16
     $               icxp,icyp,iczp, idxp,idyp,idzp, 1,ldeb,fmul,bc,ibc)11d14s22
         do i=0,ntemp-1                                                   10d8s15
          bc(itemp+i)=bc(itemp+i)+bc(itemp2+i)                          10d8s15
         end do                                                         10d8s15
        end do                                                          10d8s15
       end do                                                           10d8s15
       jtemp=itemp                                                      10d8s15
       do ia=0,nsza-1                                                   10d8s15
        iap=ia+nsza
        do ib=0,nszb-1                                                  10d8s15
         ibp=ib+nsza
         do ic=0,nszc-1                                                 10d8s15
          icp=ic+nszc                                                   10d8s15
          do id=0,nszd-1                                                10d8s15
           jd=id+icartt                                                 10d8s15
           jdp=jd+nszd                                                  10d8s15
           iada=jdp+nszdc*(icp+nszcc*(ibp+nszbc*iap))                     10d8s15
           bc(iada)=bc(jtemp)                                           10d8s15
           jtemp=jtemp+1                                                10d8s15
          end do                                                        10d8s15
         end do                                                         10d8s15
        end do                                                          10d8s15
       end do                                                           10d8s15
       ibcoff=itemp                                                     10d8s15
      end if
      if(idorel.lt.0)then
c
c     breit correction
c
       ntemp=nsza*nszb*nszc*nszd                                        10d8s15
       itemp=ibcoff                                                     10d8s15
       ibcoff=itemp+ntemp                                               10d9s15
       call enough('restid.  4',bc,ibc)
       do k=1,4
        do i=0,ntemp-1                                                   10d8s15
         bc(itemp+i)=0d0                                                10d8s15
        end do                                                           10d8s15
        do ib=1,12                                                       10d8s15
         fmul=ascale*breit(ib)
         ldeb=.false.
         iaxp=iax+ibt(1,ib,k)                                           4d27s16
         iayp=iay+ibt(2,ib,k)                                           4d27s16
         iazp=iaz+ibt(3,ib,k)                                           4d27s16
         ibxp=ibx+ibt(4,ib,k)                                           4d27s16
         ibyp=iby+ibt(5,ib,k)                                           4d27s16
         ibzp=ibz+ibt(6,ib,k)                                           4d27s16
         icxp=icx+ibt(7,ib,k)                                           4d27s16
         icyp=icy+ibt(8,ib,k)                                           4d27s16
         iczp=icz+ibt(9,ib,k)                                           4d27s16
         idxp=idx+ibt(10,ib,k)                                           4d27s16
         idyp=idy+ibt(11,ib,k)                                           4d27s16
         idzp=idz+ibt(12,ib,k)                                           4d27s16
         call derid(la,zza,xna,xa,ya,za,iadda,                            10d8s15
     $               lb,zzb,xnb,xb,yb,zb,iaddb,                         1d29s10
     $               lc,zzc,xnc,xc,yc,zc,iaddc,                         1d29s10
     $               ld,zzd,xnd,xd,yd,zd,iaddd,eraw,nbasall,nsymb,      5d11s10
     $               itemp2,iaxp,iayp,iazp, ibxp,ibyp,ibzp,             4d27s16
     $               icxp,icyp,iczp, idxp,idyp,idzp,ibt(13,ib,k),ldeb,  4d27s16
     $        fmul,bc,ibc)                                              11d14s22
         do i=0,ntemp-1                                                   10d8s15
          bc(itemp+i)=bc(itemp+i)+bc(itemp2+i)                          10d9s15
         end do                                                         10d8s15
        end do                                                          10d8s15
        nap=nsza*(ibt(1,1,k)+ibt(2,1,k)+ibt(3,1,k))
        nbp=nszb*(ibt(4,1,k)+ibt(5,1,k)+ibt(6,1,k))
        ncp=nszc*(ibt(7,1,k)+ibt(8,1,k)+ibt(9,1,k))
        ndp=nszd*(ibt(10,1,k)+ibt(11,1,k)+ibt(12,1,k))
        jtemp=itemp                                                     10d9s15
        do ia=0,nsza-1                                                   10d8s15
         iap=ia+nap                                                     10d9s15
         do ib=0,nszb-1                                                  10d8s15
          ibp=ib+nbp                                                    10d9s15
          do ic=0,nszc-1                                                 10d8s15
           icp=ic+ncp                                                   10d9s15
           do id=0,nszd-1                                                10d8s15
            jdp=id+icartt+ndp                                           10d9s15
            iada=jdp+nszdc*(icp+nszcc*(ibp+nszbc*iap))                  10d9s15
            bc(iada)=bc(jtemp)                                           10d8s15
            jtemp=jtemp+1                                                10d8s15
           end do                                                        10d8s15
          end do                                                         10d8s15
         end do                                                          10d8s15
        end do                                                           10d8s15
       end do                                                           10d9s15
       ibcoff=itemp                                                     10d8s15
      end if                                                            10d8s15
      ibcoff=ibcoffo                                                    10d8s15
      return
      end
