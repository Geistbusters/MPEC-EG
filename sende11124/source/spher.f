c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine spher(lmax,bc,ibc)                                     11d9s22
      implicit real*8 (a-h,o-z)
c
c     compute cartesian to spherical transformations
c
      parameter (id=40)
      parameter (id2=id*2)
      character*3 stype(8)
      dimension ctheta(id),wtheta(id),phi(id2),wphi(id2),tmp(id),
     $     stheta(id),ipow(4,id2),fact(id,3),fcn(id2),plm(id),
     $     sph(id,id2),cph(id,id2),coef(id,id2),xnorm(id),
     $     plmu(id),isym(3,8),nhere(8),mhere(8),ml(id,8),
     $     iptl(id,8),smat(id2,id2),coeft(id2,id),ipvt(id2),
     $     nnz(id),iptnz(id),jptnz(id),kptnz(id)                        11d3s16
      include "common.spher"                                            2d19s10
      include "common.store"                                            2d19s10
      data isym/0,0,0, 1,0,0, 0,1,0, 1,1,0, 0,0,1, 1,0,1, 0,1,1, 1,1,1/
      data stype/'Ag ','B3u','B2u','B1g','B1u','B2g','B3g','Au '/
      pi=acos(-1d0)
      fact(1,1)=1d0
      fact(1,2)=1d0
      fact(1,3)=1d0
      do l=0,lmax
       lp=l+1                                                           2d19s10
       do i=1,8                                                         2d19s10
        ipt(i,lp)=0                                                     2d19s10
       end do                                                           2d19s10
       ip=0
       do ipass=1,8
        iwrt=0
        nhere(ipass)=0
        mhere(ipass)=0
        do ix=0,l
         iyx=l-ix
         do iy=0,iyx
          iz=l-ix-iy
          if(mod(ix,2).eq.isym(1,ipass).and.mod(iy,2).eq.isym(2,ipass)
     $         .and.mod(iz,2).eq.isym(3,ipass))then
           if(iwrt.eq.0)then
    4       format('symmetry: ',a3)
            iwrt=1
           end if
           ip=ip+1
           if(ip.gt.id2)stop 'ip'
           ipow(1,ip)=ix+1
           ipow(2,ip)=iy+1
           ipow(3,ip)=iz+1
           ipow(4,ip)=ipass
           nhere(ipass)=nhere(ipass)+1
    1      format(i5,5x,3i3)
          end if
         end do
        end do
       end do
       ipp(lp)=ibcoff                                                   2d22s10
       ibcoff=ipp(lp)+3*ip+1                                            2d22s10
       call enough('spher.  1',bc,ibc)
       jj=ipp(lp)+1                                                     2d22s10
       do i=1,ip                                                        2d22s10
        ibc(jj)=ipow(1,i)-1                                             2d22s10
        jj=jj+1                                                         2d22s10
        ibc(jj)=ipow(2,i)-1                                             2d22s10
        jj=jj+1                                                         2d22s10
        ibc(jj)=ipow(3,i)-1                                             2d22s10
        jj=jj+1                                                         2d22s10
       end do                                                           2d22s10
       ibc(ipp(lp))=ip                                                  2d22s10
       num=2*l+1
       if(num.gt.id)stop 'num'
       do i=1,ip
        do j=1,num
         coef(j,i)=0d0
        end do
        do j=1,ip
         smat(j,i)=0d0
        end do
       end do
       do i=1,num
        xnorm(i)=0d0
       end do
       nq=l+2
       call gaussq(1,nq,dum,dum,0,dum,tmp,ctheta,wtheta)
       call gaussq(2,nq,dum,dum,0,dum,tmp,phi,wphi)
       nq2=nq*2
       sum1=0d0
       sum2=0d0
       do i=1,nq
        sum1=sum1+wtheta(i)
        sum2=sum2+wphi(i)
        stheta(i)=sqrt(1d0-ctheta(i)**2)
        phi(i)=acos(phi(i))
        phi(i+nq)=phi(i)+pi
        wphi(i+nq)=wphi(i)
       end do
       do i=1,nq2
        do m=0,max(1,l)
         sph(m+1,i)=sin(dfloat(m)*phi(i))
         cph(m+1,i)=cos(dfloat(m)*phi(i))
        end do
       end do
       do it=1,nq
        do m=0,l
         call plmv(tmp,l,m,ctheta(it),1)
         plm(m+1)=tmp(l+1-m)*wtheta(it)
         plmu(m+1)=tmp(l+1-m)
        end do
        do iph=1,nq2
         sp=sph(2,iph)
         cp=cph(2,iph)
         x=stheta(it)*cp
         y=stheta(it)*sp
         z=ctheta(it)
         fact(1,1)=wphi(iph)
         do i=1,l
          fact(i+1,1)=fact(i,1)*x
          fact(i+1,2)=fact(i,2)*y
          fact(i+1,3)=fact(i,3)*z
         end do
         do ifcn=1,ip
          fcn(ifcn)=fact(ipow(1,ifcn),1)*fact(ipow(2,ifcn),2)
     $         *fact(ipow(3,ifcn),3)
         end do
         do ifcn=1,ip
          xx=fcn(ifcn)*wtheta(it)/wphi(iph)
          do jfcn=1,ip
           smat(jfcn,ifcn)=smat(jfcn,ifcn)+fcn(jfcn)*xx
          end do
         end do
         xnorm(1)=xnorm(1)+plmu(1)*wphi(iph)*plm(1)
          mm=2
          do m=1,l
           mp=m+1
           term=plmu(mp)*plm(mp)*wphi(iph)
           xnorm(mm)=xnorm(mm)+term*(cph(mp,iph)**2)
           mm=mm+1
           xnorm(mm)=xnorm(mm)+term*(sph(mp,iph)**2)
           mm=mm+1
          end do
          do ifcn=1,ip
           coef(1,ifcn)=coef(1,ifcn)+fcn(ifcn)*plm(1)
           mm=2
           do m=1,l
            mp=m+1
            term=fcn(ifcn)*plm(mp)
            coef(mm,ifcn)=coef(mm,ifcn)+term*cph(mp,iph)
            mm=mm+1
            coef(mm,ifcn)=coef(mm,ifcn)+term*sph(mp,iph)
            mm=mm+1
           end do
          end do
        end do
       end do
       do i=1,num
        xnorm(i)=1d0/sqrt(xnorm(i))
       end do
       ihit=0
       do ifcn=1,ip
        coef(1,ifcn)=coef(1,ifcn)*xnorm(1)
        if(abs(coef(1,ifcn)).gt.1d-12)then
         if(ihit.eq.0)then
          ihit=1
          mhere(ipow(4,ifcn))=1
          ml(1,ipow(4,ifcn))=0
          iptl(1,ipow(4,ifcn))=1
         end if
        end if
       end do
    2  format('m =',i3,1pe15.7)
       mm=2
       do i=1,l
        ihit=0
        do ifcn=1,ip
         coef(mm,ifcn)=coef(mm,ifcn)*xnorm(mm)
         if(abs(coef(mm,ifcn)).gt.1d-12)then
          if(ihit.eq.0)then
           ihit=1
           mhere(ipow(4,ifcn))=mhere(ipow(4,ifcn))+1
           ml(mhere(ipow(4,ifcn)),ipow(4,ifcn))=i
           iptl(mhere(ipow(4,ifcn)),ipow(4,ifcn))=mm
          end if
         end if
        end do
        mm=mm+1
        ihit=0
        do ifcn=1,ip
         coef(mm,ifcn)=coef(mm,ifcn)*xnorm(mm)
         if(abs(coef(mm,ifcn)).gt.1d-8)then
          if(ihit.eq.0)then
           ihit=1
           mhere(ipow(4,ifcn))=mhere(ipow(4,ifcn))+1
           ml(mhere(ipow(4,ifcn)),ipow(4,ifcn))=-i
           iptl(mhere(ipow(4,ifcn)),ipow(4,ifcn))=mm
          end if
         end if
        end do
        mm=mm+1
       end do
       do i=1,ip
        do j=1,num
         coeft(i,j)=coef(j,i)
        end do
       end do
       call lusolv(smat,id2,ip,coeft,id2,num,ipvt,ierr,3)
       if(ierr.ne.0)then
        write(6,*)('in sphersum, on return from lusolv, ierr = '),ierr
        stop
       end if
       do i=1,ip
        do j=1,num
         coef(j,i)=coeft(i,j)
        end do
       end do
   32  format('for l =',i3)
       nnzx=0                                                           11d3s16
       do i=1,id                                                        11d3s16
        nnz(i)=0                                                        11d3s16
       end do                                                           11d3s16
       npowtot=0                                                        2d19s10
       jmlzb=0                                                          4d19s21
       do jsym=1,8
        if(nhere(jsym).gt.0)then
         npowtot=npowtot+nhere(jsym)                                    2d19s10
         ipt(jsym,lp)=ibcoff                                            2d19s10
         ibcoff=ibcoff+2+mhere(jsym)*(nhere(jsym)+1)+nhere(jsym)*3      2d19s10
         call enough('spher.  2',bc,ibc)
         ibc(ipt(jsym,lp))=nhere(jsym)                                  2d19s10
         ibc(ipt(jsym,lp)+1)=mhere(jsym)                                2d19s10
         do i=1,mhere(jsym)                                             2d19s10
          ibc(ipt(jsym,lp)+1+i)=ml(i,jsym)                              2d19s10
          mlzb(jmlzb+i,lp)=iabs(ml(i,jsym))                             4d19s21
         end do                                                         2d19s10
         jmlzb=jmlzb+mhere(jsym)                                        4d19s21
         ipowr=ipt(jsym,lp)+2+mhere(jsym)                                2d19s10
         jpowr=ipowr                                                      2d19s10
         icoef=ipowr+nhere(jsym)*3                                       2d19s10
   30    format(i2,2i3)
   31    format(10i3)
         jcoef=icoef-nhere(jsym)                                        2d19s10
         do ifcn=1,ip
          if(ipow(4,ifcn).eq.jsym)then
   22      format(3i3,1p10e25.17)
           ibc(jpowr)=ipow(1,ifcn)-1                                    2d19s10
           ibc(jpowr+1)=ipow(2,ifcn)-1                                  2d19s10
           ibc(jpowr+2)=ipow(3,ifcn)-1                                  2d19s10
           jpowr=jpowr+3                                                  2d19s10
           do i=1,mhere(jsym)                                           2d19s10
            bc(jcoef+i*nhere(jsym))=coef(iptl(i,jsym),ifcn)             2d19s10
           end do                                                       2d19s10
           jcoef=jcoef+1                                                2d19s10
          end if
         end do
        do i=1,mhere(jsym)
         nz=0
         sz=0d0
         do j=1,nhere(jsym)
          ji=icoef+j-1+nhere(jsym)*(i-1)
          sz=sz+bc(ji)**2
         end do
         sz=sqrt(sz/dfloat(nhere(jsym)))
         nx=0
         do j=1,nhere(jsym)
          ji=icoef+j-1+nhere(jsym)*(i-1)
          if(abs(bc(ji)).gt.sz*1d-10)nz=nz+1
         end do
         nnzx=max(nnzx,nz)                                              11d3s16
         if(nz.le.id)then                                               11d3s16
          nnz(nz)=nnz(nz)+1                                             11d3s16
         end if                                                         11d3s16
        end do
   20   format(20i5)
        end if
       end do
       ioff=1                                                           11d3s16
       ncf=0                                                            11d3s16
       ndo=0                                                            11d3s16
       do i=1,min(nnzx,id)                                              11d3s16
        ncf=ncf+nnz(i)*i                                                11d4s16
        iptnz(i)=ioff                                                   11d3s16
        if(nnz(i).gt.0)then                                             11d3s16
         ndo=ndo+1                                                      11d3s16
         nstor=3+nnz(i)*(i+1)                                           11d4s16
         ioff=ioff+nstor                                                11d3s16
        end if
       end do                                                           11d3s16
       idat1=ibcoff                                                     11d3s16
       idat2=idat1+ioff                                                 11d3s16
       ibcoff=idat2+ncf                                                 11d3s16
       call enough('spher.  3',bc,ibc)
       idat1=idat1*intmul                                               11d3s16
       jdat2=idat2                                                      11d3s16
       ibc(idat1)=ndo                                                   11d3s16
       do i=1,min(nnzx,id)                                              11d3s16
        iptnz(i)=iptnz(i)+idat1                                         11d3s16
        jptnz(i)=iptnz(i)                                               11d3s16
        if(nnz(i).gt.0)then                                             11d3s16
         ibc(jptnz(i))=i                                                11d3s16
         jptnz(i)=jptnz(i)+1                                            11d3s16
         ibc(jptnz(i))=nnz(i)                                           11d3s16
         jptnz(i)=jptnz(i)+1                                            11d3s16
         ibc(jptnz(i))=jdat2                                            11d3s16
         jptnz(i)=jptnz(i)+1                                            11d3s16
        end if                                                          11d3s16
        kptnz(i)=jdat2                                                  11d3s16
        jdat2=jdat2+nnz(i)*i                                            11d3s16
       end do                                                           11d3s16
       noff=0                                                           11d3s16
       moff=0                                                           11d3s16
       do jsym=1,8                                                      11d3s16
        if(nhere(jsym).gt.0)then                                        11d3s16
         ipowr=ipt(jsym,lp)+2+mhere(jsym)
         icoef=ipowr+nhere(jsym)*3
         do i=1,mhere(jsym)                                             11d3s16
          sz=0d0
          do j=1,nhere(jsym)
           ji=icoef+j-1+nhere(jsym)*(i-1)
           sz=sz+bc(ji)**2
          end do
          sz=sqrt(sz/dfloat(nhere(jsym)))
          nz=0
          do j=1,nhere(jsym)
           ji=icoef+j-1+nhere(jsym)*(i-1)
           if(abs(bc(ji)).gt.sz*1d-10)nz=nz+1
          end do
          do j=1,nhere(jsym)
           ji=icoef+j-1+nhere(jsym)*(i-1)
           if(abs(bc(ji)).gt.sz*1d-10)then
            bc(kptnz(nz))=bc(ji)
            kptnz(nz)=kptnz(nz)+1
            ibc(jptnz(nz))=noff+j                                       11d3s16
            jptnz(nz)=jptnz(nz)+1                                       11d3s16
           end if
          end do
          ibc(jptnz(nz))=i+moff                                         11d3s16
          jptnz(nz)=jptnz(nz)+1                                         11d3s16
         end do                                                         11d3s16
         noff=noff+nhere(jsym)                                          11d3s16
         moff=moff+mhere(jsym)                                          11d3s16
        end if                                                          11d3s16
       end do                                                           11d3s16
       jj=idat1+1
       ippn(lp)=idat1                                                   11d4s16
       do i=1,ibc(idat1)
        nt=ibc(jj)
        jj=jj+1
        nz=ibc(jj)
        jj=jj+1
        kk=ibc(jj)
        jj=jj+1
        do j=1,nz
         ito=ibc(jj+nt)
         do k=1,nt
          jj=jj+1
          kk=kk+1
         end do
         jj=jj+1
        end do
       end do
   33  format('eol')
      end do
      return
      end
