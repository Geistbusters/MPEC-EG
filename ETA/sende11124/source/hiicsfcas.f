c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hiicsfcas(hdiag,nfcn,mdon,mdoo,ibasis,                 8d2s22
     $     ih0a,i2e,shift,nec,ism,irel,irefo,iptrbit,mysym,norb,bc,ibc) 11d10s22
c
c     compute spin independent part of hii
c
      implicit real*8 (a-h,o-z)
      logical ldebug                                                    7d1s19
      dimension hdiag(nfcn),ih0a(*),i2e(*),ibasis(*),ism(*),            8d2s22
     $     irel(*),irefo(*),iptrbit(2,mdoo+1,*)                         8d2s22
      include "common.store"
      include "common.print"                                            1d23s19
      if(iprtr(12).eq.0)then                                            1d23s19
       ldebug=.false.                                                    7d15s19
      else                                                              1d23s19
       ldebug=.true.
      end if                                                            1d23s19
      if(ldebug)write(6,*)('hi, i''m hiicsfcas'),nfcn                           7d1s19
      if(ldebug)then
       write(6,*)('h0: ')                                               1d24s20
       do i=1,8                                                         1d24s20
        if(irefo(i).gt.0)then                                           1d24s20
         write(6,*)('for symmetry block '),i,ih0a(i)                            1d24s20
         call prntm2(bc(ih0a(i)),irefo(i),irefo(i),irefo(i))            1d24s20
        end if                                                          1d24s20
       end do
      end if                                                            7d15s19
      nclo=-1
      ihdiag=1
      do if=1,nfcn
       jbasis=3*(if-1)
       if(ldebug)                                                       7d1s19
     $    write(6,*)('fcn no. '),if,('basis: '),(ibasis(jbasis+j),j=1,3)7d1s19
       iic=iptrbit(1,ibasis(jbasis+1)+1,mysym)+ibasis(jbasis+2)-1       8d2s22
       iio=iptrbit(2,ibasis(jbasis+1)+1,mysym)+ibasis(jbasis+3)-1       8d2s22
       if(ldebug)then
        call dcbit(ibc(iic),norb,'iic')
        call dcbit(ibc(iio),norb,'iio')
       end if
       xc=0d0                                                           8d2s22
       xo=0d0                                                           8d2s22
       do i=1,norb                                                      8d2s22
        isb=ism(i)                                                      8d2s22
        ig=irel(i)-1                                                    8d2s22
        ih=ih0a(isb)+ig*(irefo(isb)+1)                                  8d2s22
        njrow=(irefo(isb)*(irefo(isb)+1))/2                             8d2s22
        ijrow=ig+((ig*(ig+1))/2)
        if(btest(ibc(iic),i))then                                       8d2s22
         orig=xc
         xc=xc+2d0*bc(ih)                                               8d2s22
         do j=1,norb                                                    8d2s22
          if(btest(ibc(iic),j))then                                     8d2s22
           jsb=ism(j)                                                   8d2s22
           jg=irel(j)-1                                                 8d2s22
           i2eu=ifind2(isb,isb,jsb,jsb,icase)
           ijcol=jg+((jg*(jg+1))/2)
           iad=i2e(i2eu)+ijrow+njrow*ijcol
           orig=xc
           xc=xc+2d0*bc(iad)                                            8d2s22
           i2eu=ifind2(isb,jsb,isb,jsb,icase)                           6d14s19
c                                                                       6d14s19
c     icase=1 1234   ig jg ig jg                                        6d14s19
c     icase=2 2143   jg ig jg ig                                        6d14s19
c     icase=3 1243   ig jg jg ig                                        6d14s19
c     icase=4 2134   jg ig ig jg                                        6d14s19
c                                                                       6d14s19
           if(isb.eq.jsb)then                                           6d14s19
            ix=max(jg,ig)                                               6d14s19
            in=min(jg,ig)                                               6d14s19
            iad=i2e(i2eu)+(((ix*(ix+1))/2)+in)*(njrow+1)                6d14s19
           else if(icase.eq.1)then                                      6d14s19
            iad=i2e(i2eu)+ig+irefo(isb)*(jg+irefo(jsb)                  6d14s19
     $           *(ig+irefo(isb)*jg))                                   6d14s19
           else if(icase.eq.2)then                                      6d14s19
            iad=i2e(i2eu)+jg+irefo(jsb)*(ig+irefo(isb)                  6d14s19
     $           *(jg+irefo(jsb)*ig))                                   6d14s19
           else if(icase.eq.3)then                                      6d14s19
            iad=i2e(i2eu)+ig+irefo(isb)*(jg+irefo(jsb)                  6d14s19
     $           *(jg+irefo(jsb)*ig))                                   6d14s19
           else                                                         6d14s19
            iad=i2e(i2eu)+jg+irefo(jsb)*(ig+irefo(isb)                  6d14s19
     $           *(ig+irefo(isb)*jg))                                   6d14s19
           end if                                                       6d14s19
           orig=xc
           xc=xc-bc(iad)                                                8d2s22
          end if                                                        8d2s22
         end do                                                         8d2s22
        end if                                                          8d2s22
        if(btest(ibc(iio),i))then                                       8d2s22
         xo=xo+bc(ih)                                                   8d2s22
         do j=1,norb                                                    8d2s22
          if(btest(ibc(iio),j).and.j.ne.i)then                          8d2s22
           jsb=ism(j)                                                   8d2s22
           jg=irel(j)-1                                                 8d2s22
           i2eu=ifind2(isb,isb,jsb,jsb,icase)                           8d2s22
           ijcol=jg+((jg*(jg+1))/2)                                     8d2s22
           iad=i2e(i2eu)+ijrow+njrow*ijcol                              8d2s22
           xo=xo+0.5d0*bc(iad)                                          8d2s22
          end if                                                        8d2s22
         end do                                                         8d2s22
        end if                                                          8d2s22
       end do                                                           8d2s22
       nclo=ibasis(jbasis+1)                                            8d2s22
       nclop=nclo+1                                                     8d2s22
       nopen=nec-2*nclo                                                 8d2s22
       xd=xc+xo                                                         8d2s22
       hd=xd                                                            8d2s22
       if(ldebug)write(6,*)('sum of oo and cc parts: '),hd,xc,xo              7d1s19
       hd=hd+shift
       if(ldebug)write(6,*)('with shift '),hd                           7d1s19
       jjc=iptrbit(1,nclop,mysym)+ibasis(jbasis+2)-1                    8d2s22
       jjo=iptrbit(2,nclop,mysym)+ibasis(jbasis+3)-1                    8d2s22
       do i=1,norb                                                      8d2s22
        if(btest(ibc(jjc),i))then                                       8d2s22
         iorb=i                                                         8d2s22
         isb=ism(iorb)
         ig=irel(iorb)-1                                                 6d11s19
         njrow=(irefo(isb)*(irefo(isb)+1))/2                             6d11s19
         ijrow=ig+((ig*(ig+1))/2)                                        6d11s19
         do j=1,norb                                                     8d2s22
          if(btest(ibc(jjo),j))then                                      8d2s22
           jorb=j                                                        8d2s22
           jsb=ism(jorb)                                                  6d11s19
           jg=irel(jorb)-1                                                6d11s19
           ijcol=jg+((jg*(jg+1))/2)                                       6d11s19
           i2eu=ifind2(isb,isb,jsb,jsb,icase)                             6d11s19
           iad=i2e(i2eu)+ijrow+njrow*ijcol                                6d11s19
           hd=hd+2d0*bc(iad)                                              6d11s19
           if(ldebug)write(6,*)('J: '),bc(iad),i2eu,isb,jsb
           i2eu=ifind2(isb,jsb,isb,jsb,icase)                             6d11s19
           if(isb.eq.jsb)then                                             6d11s19
            ix=max(jg,ig)                                                 6d11s19
            in=min(jg,ig)                                                 6d11s19
            iad=i2e(i2eu)+(in+((ix*(ix+1))/2))*(njrow+1)                  6d14s19
           else if(icase.eq.1)then                                      6d14s19
            iad=i2e(i2eu)+ig+irefo(isb)*(jg+irefo(jsb)                  6d14s19
     $           *(ig+irefo(isb)*jg))                                   6d14s19
           else if(icase.eq.2)then                                      6d14s19
            iad=i2e(i2eu)+jg+irefo(jsb)*(ig+irefo(isb)                  6d14s19
     $           *(jg+irefo(jsb)*ig))                                   6d14s19
           else if(icase.eq.3)then                                      6d14s19
            iad=i2e(i2eu)+ig+irefo(isb)*(jg+irefo(jsb)                  6d14s19
     $           *(jg+irefo(jsb)*ig))                                   6d14s19
           else                                                         6d14s19
            iad=i2e(i2eu)+jg+irefo(jsb)*(ig+irefo(isb)                  6d14s19
     $           *(ig+irefo(isb)*jg))                                   6d14s19
           end if                                                       6d14s19
           if(ldebug)write(6,*)('K: '),bc(iad)                            1d24s20
           hd=hd-bc(iad)
          end if                                                         8d2s22
         end do                                                          6d11s19
        end if                                                          8d2s22
       end do
       if(ldebug)write(6,*)('final spin independent part: '),hd         7d1s19
       hdiag(ihdiag)=hd
       ihdiag=ihdiag+1
      end do
      return
      end
