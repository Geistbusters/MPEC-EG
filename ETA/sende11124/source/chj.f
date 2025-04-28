c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine chj(h,n,eig,vec,eposcut,nes,iden,nvec,unvecs,nstate,   2d26s20
     $     iocc,iflag,scopy,nbaslarge,nbassmall,bc,ibc)                 11d9s22
      implicit real*8 (a-h,o-y)                                         2d22s20
      implicit complex*16 (z)
      logical lprint
c
c     diagonalize complex hermitian matrix h by a complex version of    2d22s20
c     the jacobi method.
c     iflag=0 means just diagonalize
c
      integer*8 iocc(nvec,nstate)                                       2d26s20
      dimension h(n,n,2),eig(n),vec(n,n,2),zv(2,2),zh(2,2),zt(2,2),
     $     unvecs(n,n),scopy(*)                                         3d16s20
      include "common.store"
      include "common.basis"
      include "common.input"
      ibcoffo=ibcoff
      lprint=.false.
      n2=n*n
      thrs=-smallest*eposcut*dfloat(n)                                  2d25s20
      srh=sqrt(0.5d0)                                                   2d22s20
      pi=acos(-1d0)
      szr=0d0
      szi=0d0
      do i=1,n-1
       do j=i+1,n
        diffr=abs(h(j,i,1)-h(i,j,1))
        diffi=abs(h(j,i,2)+h(i,j,2))
        szr=szr+(h(j,i,1)-h(i,j,1))**2
        szi=szi+(h(j,i,2)+h(i,j,2))**2
       end do
      end do
      szr=sqrt(szr/dfloat(n*n))
      szi=sqrt(szi/dfloat(n*n))
      if(max(szr,szi).gt.thrs)then
       write(6,*)('hermitian test results: '),szr,szi
       write(6,*)('vs. thrs = '),thrs
       call dws_sync
       call dws_finalize
       stop 'chj'
      end if
      ihcpy=ibcoff
      ibcoff=ihcpy+n*n*2
      call enough('chj.  1',bc,ibc)
      jhcpy=ihcpy-1
      do i3=1,2                                                         1d17s23
       do i2=1,n                                                        1d17s23
        do i1=1,n                                                       1d17s23
         bc(jhcpy+i1)=h(i1,i2,i3)                                       1d17s23
        end do                                                          1d17s23
        jhcpy=jhcpy+n                                                   1d17s23
       end do                                                           1d17s23
      end do                                                            1d17s23
      do i=1,n                                                          2d22s20
       do j=1,n                                                         2d22s20
        vec(j,i,1)=0d0                                                  2d22s20
        vec(j,i,2)=0d0                                                  2d22s20
       end do                                                           2d22s20
       vec(i,i,1)=1d0                                                   2d22s20
      end do                                                            2d22s20
      isweep=0                                                          2d22s20
    1 continue                                                          2d22s20
       isweep=isweep+1                                                  2d22s20
       if(isweep.gt.1000)then
        write(6,*)('real part so far: ')
        call prntm2(h,n,n,n)
        write(6,*)('imag part so far: ')
        call prntm2(h(1,1,2),n,n,n)
      szr=0d0
      szi=0d0
      do i=1,n-1
       do j=i+1,n
        diffr=abs(h(j,i,1)-h(i,j,1))
        diffi=abs(h(j,i,2)+h(i,j,2))
        if(max(diffr,diffi).gt.thrs)write(6,*)j,i,h(j,i,1),h(i,j,1),
     $       h(j,i,2),h(i,j,2)
        szr=szr+(h(j,i,1)-h(i,j,1))**2
        szi=szi+(h(j,i,2)+h(i,j,2))**2
       end do
      end do
      szr=sqrt(szr/dfloat(n*n))
      szi=sqrt(szi/dfloat(n*n))
      write(6,*)('sym test: '),szr,szi
      write(6,*)('transform original by vectors: ')
      write(6,*)('original real: ')
      call prntm2(bc(ihcpy),n,n,n)
      write(6,*)('imag')
      call prntm2(bc(ihcpy+n2),n,n,n)
      write(6,*)('real part of vectors: ')
      call prntm2(vec,n,n,n)
      write(6,*)('imag ')
      call prntm2(vec(1,1,2),n,n,n)
      itmp=ibcoff
      ibcoff=itmp+n*n*2
      call enough('chj.  2',bc,ibc)
      write(6,*)('unitarity test: ')
      call dgemm('t','n',n,n,n,1d0,vec,n,vec,n,0d0,bc(itmp),n,
     d' chj.  1')
      call dgemm('t','n',n,n,n,1d0,vec(1,1,2),n,vec(1,1,2),n,1d0,
     $     bc(itmp),n,
     d' chj.  2')
      itmi=itmp+n2
      call dgemm('t','n',n,n,n,1d0,vec,n,vec(1,1,2),n,0d0,bc(itmi),n,
     d' chj.  3')
      call dgemm('t','n',n,n,n,-1d0,vec(1,1,2),n,vec,n,1d0,bc(itmi),n,
     d' chj.  4')
      write(6,*)('real part: ')
      call prntm2(bc(itmp),n,n,n)
      write(6,*)('imag part: ')
      call prntm2(bc(itmi),n,n,n)
      ihcpi=ihcpy+n2
      do ipass=1,2
       call dgemm('n','n',n,n,n,1d0,bc(ihcpy),n,vec,n,0d0,bc(itmp),n,
     d' chj.  5')
       call dgemm('n','n',n,n,n,-1d0,bc(ihcpi),n,vec(1,1,2),n,1d0,
     $      bc(itmp),n,
     d' chj.  6')
       call dgemm('n','n',n,n,n,1d0,bc(ihcpy),n,vec(1,1,2),n,0d0,
     $      bc(itmi),n,
     d' chj.  7')
       call dgemm('n','n',n,n,n,1d0,bc(ihcpi),n,vec,n,1d0,bc(itmi),n,
     d' chj.  8')
       do i=0,n-1
        do j=0,n-1
         jir=itmp+j+n*i
         jii=jir+n2
         ijr=ihcpy+i+n*j
         iji=ijr+n2
         bc(ijr)=bc(jir)
         bc(iji)=-bc(jii)
        end do
       end do
      end do
      write(6,*)('real part: ')
      call prntm2(bc(ihcpy),n,n,n)
      write(6,*)('imag part: ')
      call prntm2(bc(ihcpy+n2),n,n,n)
        call dws_sync
        call dws_finalize
        stop 'isweep'
       end if
       hx=0d0                                                           2d22s20
       do i=1,n-1                                                       2d22s20
        do j=i+1,n                                                      2d22s20
         hx=max(abs(h(j,i,1)),abs(h(j,i,2)),hx)                         2d22s20
        end do                                                          2d22s20
       end do                                                           2d22s20
       if(hx.gt.thrs)then
        thrsh=hx*0.1d0
        do i=1,n-1                                                      2d22s20
         do j=i+1,n                                                     2d22s20
          if(lprint)then
           write(6,*)('for j,i '),j,i
           write(6,*)('diags: '),h(i,i,1),h(j,j,1)
           write(6,*)('off '),h(j,i,1),h(j,i,2)
          end if
          sz=sqrt(h(j,i,1)**2+h(j,i,2)**2)
          if(sz.gt.thrsh)then
           b=h(j,j,1)-h(i,i,1)
           if(abs(b).lt.thrs)then                                       2d23s20
            a=srh                                                       2d23s20
            br=h(j,i,1)**2-h(j,i,2)**2                                  2d23s20
            bi=2d0*h(j,i,1)*h(j,i,2)                                    2d22s20
            sz=1d0/sqrt(br**2+bi**2)                                    2d23s20
            br=srh*br*sz                                                2d23s20
            bi=srh*bi*sz                                                2d23s20
           else                                                         2d23s20
            acm=h(j,i,1)**2+h(j,i,2)**2                                   2d22s20
            desc=sqrt(1d0+4d0*(acm/(b*b)))                                2d22s20
            if(lprint)write(6,*)('descriminant: '),desc                             2d22s20
            dp=1d0+desc                                                   2d22s20
            dm=1d0-desc                                                   2d22s20
            xfactr=-0.5d0*b/(h(j,i,1)**2+h(j,i,2)**2)                     2d22s20
            xfacti=xfactr*h(j,i,2)                                        2d22s20
            xfactr=xfactr*h(j,i,1)                                        2d22s20
            x1r=dp*xfactr
            x1i=dp*xfacti
            x2r=dm*xfactr
            x2i=dm*xfacti
            if(lprint)write(6,*)('x roots: '),x1r,x1i,x2r,x2i
            xi1r=-dp*xfactr
            xi1i=dp*xfacti
            xi2r=-dm*xfactr
            xi2i=dm*xfacti
            if(lprint)write(6,*)('xi roots: '),xi1r,xi1i,xi2r,xi2i
            xii1r=xi1r/(xi1r**2+xi1i**2)                                  2d22s20
            xii1i=-xi1i/(xi1r**2+xi1i**2)                                  2d22s20
            xii2r=xi2r/(max(1d-10,xi2r**2+xi2i**2))                     1d17s23
            xii2i=-xi2i/(max(1d-10,xi2r**2+xi2i**2))                    1d17s23
            if(lprint)write(6,*)('inverted: '),xii1r,xii1i,xii2r,xii2i
            sz1=sqrt(xii1r**2+xii1i**2)
            sz2=sqrt(x2r**2+x2i**2)
            if(lprint)write(6,*)('sizes: '),sz1,sz2
            if(sz1.gt.sz2)then                                            2d22s20
             xr=xii1r
             xi=xii1i
            else
             xr=x2r
             xi=x2i
            end if
            ts=max(sz1,sz2)
            ts2=ts*ts
            a=sqrt(ts2/(1d0+ts2))
            if(lprint)write(6,*)('a: '),a
            if(a.ne.a)then
             write(6,*)('got nan ')
             call dws_sync
             call dws_finalize
             stop
            end if
            br=a*xr/ts2
            bi=+a*xi/ts2
           end if                                                       2d23s20
           if(lprint)then
            write(6,*)('b: '),br,bi
            write(6,*)('ortho test: '),a*a+br**2+bi**2-1d0
           end if
           if(br.ne.br.or.bi.ne.bi)then
            write(6,*)('got nan ')
            call dws_sync
            call dws_finalize
            stop
           end if
           if(lprint)then
            zv(1,1)=dcmplx(a,0d0)
            zv(2,2)=zv(1,1)
            zv(2,1)=dcmplx(br,bi)
            zv(1,2)=dcmplx(-br,bi)
            write(6,*)('transformation matrix: ')
            write(6,*)(zv(1,k),k=1,2)
            write(6,*)(zv(2,k),k=1,2)
            zh(1,1)=dcmplx(h(i,i,1),0d0)
            zh(2,1)=dcmplx(h(j,i,1),h(j,i,2))
            zh(1,2)=dcmplx(h(i,j,1),h(i,j,2))
            zh(2,2)=dcmplx(h(j,j,1),0d0)
            write(6,*)('ham matrix: ')
            write(6,*)(zh(1,k),k=1,2)
            write(6,*)(zh(2,k),k=1,2)
            do ipass=1,2
             zt(1,1)=zh(1,1)*zv(1,1)+zh(1,2)*zv(2,1)
             zt(1,2)=zh(1,1)*zv(1,2)+zh(1,2)*zv(2,2)
             zt(2,1)=zh(2,1)*zv(1,1)+zh(2,2)*zv(2,1)
             zt(2,2)=zh(2,1)*zv(1,2)+zh(2,2)*zv(2,2)
             zh(1,1)=zt(1,1)
             zh(2,2)=zt(2,2)
             zh(1,2)=dconjg(zt(2,1))
             zh(2,1)=dconjg(zt(1,2))
            end do
            write(6,*)('after transformation ')
            write(6,*)(zh(1,k),k=1,2)
            write(6,*)(zh(2,k),k=1,2)
           end if
           do k=1,n
            xnewr=h(k,i,1)*a+h(k,j,1)*br-h(k,j,2)*bi
            xnewi=h(k,i,2)*a+h(k,j,1)*bi+h(k,j,2)*br
            ynewr=-br*h(k,i,1)-bi*h(k,i,2)+a*h(k,j,1)
            ynewi=bi*h(k,i,1)-br*h(k,i,2)+a*h(k,j,2)
            h(k,i,1)=xnewr
            h(k,i,2)=xnewi
            h(k,j,1)=ynewr
            h(k,j,2)=ynewi
            xnewr=vec(k,i,1)*a+vec(k,j,1)*br-vec(k,j,2)*bi
            xnewi=vec(k,i,2)*a+vec(k,j,1)*bi+vec(k,j,2)*br
            ynewr=-br*vec(k,i,1)-bi*vec(k,i,2)+a*vec(k,j,1)
            ynewi=bi*vec(k,i,1)-br*vec(k,i,2)+a*vec(k,j,2)
            vec(k,i,1)=xnewr
            vec(k,i,2)=xnewi
            vec(k,j,1)=ynewr
            vec(k,j,2)=ynewi
           end do
           do k=1,n
            xnewr=h(i,k,1)*a+h(j,k,1)*br+h(j,k,2)*bi
            xnewi=h(i,k,2)*a-h(j,k,1)*bi+h(j,k,2)*br
            ynewr=-br*h(i,k,1)+bi*h(i,k,2)+a*h(j,k,1)
            ynewi=-bi*h(i,k,1)-br*h(i,k,2)+a*h(j,k,2)
            h(i,k,1)=xnewr
            h(i,k,2)=xnewi
            h(j,k,1)=ynewr
            h(j,k,2)=ynewi
           end do
           if(lprint)then
            zh(1,1)=dcmplx(h(i,i,1),0d0)
            zh(2,1)=dcmplx(h(j,i,1),h(j,i,2))
            zh(1,2)=dcmplx(h(i,j,1),h(i,j,2))
            zh(2,2)=dcmplx(h(j,j,1),0d0)
            write(6,*)('ham matrix after transforming full matrix: ')
            write(6,*)(zh(1,k),k=1,2)
            write(6,*)(zh(2,k),k=1,2)
            lprint=.false.
           end if
          end if                                                        2d22s20
         end do                                                         2d22s20
        end do                                                          2d22s20
        go to 1
       end if
       do i=1,n
        eig(i)=h(i,i,1)
       end do
       isort=ibcoff
       ibcoff=isort+n
       call enough('chj.  3',bc,ibc)
       call dsortdws(eig,ibc(isort),n)                                  1d18s23
       ifirst=0                                                         2d28s20
       do i=1,n
        if(eig(i).le.eposcut*0.99d0)ifirst=i
       end do
       nes=n-ifirst
       ifirst=ifirst+1
       jsort=isort-1
       j=0
       ivcpy=ibcoff
       ivcpi=ivcpy+n*nes                                                2d24s20
       ibcoff=ivcpi+n*nes                                               2d24s20
       ns2=n*nes
       do i=ifirst,n
        inew=ivcpy+n*j-1
        inewi=ivcpi+n*j-1
        j=j+1
        eig(j)=eig(i)
        do k=1,n
         bc(inew+k)=vec(k,ibc(jsort+i),1)
         bc(inewi+k)=vec(k,ibc(jsort+i),2)
        end do
       end do
       if(iflag.eq.0)then                                               2d28s20
        do i=1,nes                                                      3d16s20
         iad1=ivcpy-1+n*(i-1)                                           3d16s20
         iad2=ivcpi-1+n*(i-1)                                           3d16s20
         do j=1,n                                                       3d16s20
          vec(j,i,1)=bc(iad1+j)                                         3d16s20
          vec(j,i,2)=bc(iad2+j)                                         3d16s20
         end do                                                         3d16s20
        end do                                                          3d16s20
        ibcoff=ibcoffo
        return                                                          2d28s20
       end if                                                           2d28s20
c
c     transform to ao basis
c
       call dgemm('n','n',n,nes,n,1d0,unvecs,n,bc(ivcpy),n,0d0,vec,n,   2d24s20
     d' chj.  9')
       call dgemm('n','n',n,nes,n,1d0,unvecs,n,bc(ivcpi),n,0d0,         2d24s20
     $      vec(1,1,2),n,                                               2d24s20
     d' chj. 10')
       do i=1,nes
        xmx=0d0                                                         3d22s23
        do j=1,n
         sz=vec(j,i,1)**2+vec(j,i,2)**2
         if(sz.gt.xmx)then                                              3d22s23
          xmx=sz                                                        3d22s23
          mx=j
         end if
        end do
        xmx=sqrt(xmx)                                                   3d22s23
        ca=vec(mx,i,1)/xmx                                              3d22s23
        sa=vec(mx,i,2)/xmx                                              3d22s23
        do j=1,n
         tryr=vec(j,i,1)*ca+vec(j,i,2)*sa
         tryi=-vec(j,i,1)*sa+vec(j,i,2)*ca
         vec(j,i,1)=tryr
         vec(j,i,2)=tryi
        end do
       end do
       call dws_bcast(vec,n*n*2)                                        2d26s20
       do i=0,n2*2-1                                                    2d24s20
        bc(iden+i)=0d0                                                  2d24s20
       end do                                                           2d24s20
       xavg=1d0/dfloat(nstate)                                          2d26s20
       do is=1,nstate                                                   2d26s20
        do jv=1,nvec                                                     2d24s20
         iv=iocc(jv,is)
         do i=1,n
          idenr=iden-1+n*(i-1)
          ideni=idenr+n2
          do j=1,n
           jir=ihcpy+j-1+n*(i-1)
           jii=jir+n*n
           partr=vec(j,iv,1)*vec(i,iv,1)+vec(j,iv,2)*vec(i,iv,2)
           parti=vec(j,iv,1)*vec(i,iv,2)-vec(j,iv,2)*vec(i,iv,1)
           bc(idenr+j)=bc(idenr+j)+partr*xavg                           2d26s20
           bc(ideni+j)=bc(ideni+j)+parti*xavg                           2d26s20
          end do
         end do
        end do
       end do                                                           2d26s20
      ibcoff=ibcoffo
       return
       end
