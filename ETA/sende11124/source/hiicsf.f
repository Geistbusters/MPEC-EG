c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hiicsf(hdiag,nfcn,idorb,isorb,mdon,mdoo,ibasis,iptr,   6d10s19
     $     ih0a,i2e,shift,nec,bc,ibc)                                   11d10s22
c
c     compute spin independent part of hii
c
      implicit real*8 (a-h,o-z)
      integer*1 idorb(*),isorb(*)
      logical ldebug                                                    7d1s19
      dimension hdiag(nfcn),iptr(4,*),ih0a(*),i2e(*),ibasis(*)          12d28s19
      data icall/0/                                                     1d21s21
      save icall                                                        1d21s21
      include "common.store"
      include "common.mrci"                                             6d10s19
      icall=icall+1                                                     1d21s21
      ldebug=.false.                                                     7d15s19
      if(ldebug)write(6,*)('hi, i''m hiicsf'),nfcn                           7d1s19
      if(ldebug)then
       write(6,*)('for icall = '),icall
      write(6,*)('iptr: '),loc(iptr)
      do j=1,4
       write(6,*)j,(iptr(j,k),k=1,3)
      end do
      end if                                                            7d15s19
      nclo=-1
      ihdiag=1
      igoalc=-3
      igoalcc=1750651
      igoaloo=1750653
      do if=1,nfcn
       jbasis=3*(if-1)
       if(ldebug)                                                       7d1s19
     $    write(6,*)('fcn no. '),if,('basis: '),(ibasis(jbasis+j),j=1,3)7d1s19
       if(ibasis(jbasis+1).ne.nclo)then
        if(nclo.ge.0)then
         ibcoff=ihc
        end if
        nclo=ibasis(jbasis+1)
        nclop=nclo+1
        if(ldebug)write(6,*)('rule sizes: '),iptr(1,nclop),iptr(3,nclop)7d1s19
        ihc=ibcoff
        iho=ihc+iptr(1,nclop)
        ibcoff=iho+iptr(3,nclop)
        if(nclo.eq.igoalc)write(6,*)('ihc,iho '),ihc,iho,igoalcc,igoaloo
        call enough('hiicsf.  1',bc,ibc)
        do i=0,iptr(1,nclop)+iptr(3,nclop)-1
         bc(ihc+i)=0d0
        end do
        if(ldebug)write(6,*)('for closed shells: '),iptr(2,nclop)       7d1s19
        jc=iptr(2,nclop)
        do i=0,iptr(1,nclop)-1
         ip=i+1
         if(ldebug)write(6,*)ip,(idorb(jc+j),j=0,nclo-1)                     9d20s20
         do j=0,nclo-1
          iorb=idorb(jc+j)
          isb=ism(iorb)
          ig=irel(iorb)-1
          ih=ih0a(isb)+ig*(irefo(isb)+1)                                6d10s19
          bc(ihc+i)=bc(ihc+i)+2d0*bc(ih)
          if(ldebug)write(6,*)('h0 '),iorb,bc(ih),2d0
          njrow=(irefo(isb)*(irefo(isb)+1))/2
          ijrow=ig+((ig*(ig+1))/2)                                      6d11s19
          do l=0,nclo-1
           jorb=idorb(jc+l)
           jsb=ism(jorb)
           jg=irel(jorb)-1
           i2eu=ifind2(isb,isb,jsb,jsb,icase)
           ijcol=jg+((jg*(jg+1))/2)
           iad=i2e(i2eu)+ijrow+njrow*ijcol
           bc(ihc+i)=bc(ihc+i)+2d0*bc(iad)
           if(ldebug)write(6,*)('j '),bc(iad),2d0,iorb,jorb
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
           if(ldebug)write(6,*)('k'),bc(iad),-1d0,iorb,jorb,isb,jsb,ig,
     $          jg,i2eu
           bc(ihc+i)=bc(ihc+i)-bc(iad)
          end do
         end do
         jc=jc+nclo
         if(ldebug)write(6,*)bc(ihc+i)                                  7d1s19
        end do
        nopen=nec-2*nclo                                                6d10s19
        if(ldebug)write(6,*)('for open shells: '),iptr(4,nclop)         7d1s19
        jo=iptr(4,nclop)
        do i=0,iptr(3,nclop)-1
         ip=i+1
         if(ldebug)write(6,*)ip,(isorb(jo+j),j=0,nopen-1)                    9d20s20
         do j=0,nopen-1
          iorb=isorb(jo+j)
          isb=ism(iorb)
          ig=irel(iorb)-1
          ih=ih0a(isb)+ig*(irefo(isb)+1)                                6d10s19
          bc(iho+i)=bc(iho+i)+bc(ih)
          if(ldebug)write(6,*)('h0 '),bc(ih),1d0,iorb
          njrow=(irefo(isb)*(irefo(isb)+1))/2                           6d11s19
          ijrow=ig+((ig*(ig+1))/2)
          do l=0,nopen-1
           if(l.ne.j)then                                               6d11s19
            jorb=isorb(jo+l)                                            6d11s19
            jsb=ism(jorb)                                               6d11s19
            jg=irel(jorb)-1                                             6d11s19
            i2eu=ifind2(isb,isb,jsb,jsb,icase)
            ijcol=jg+((jg*(jg+1))/2)
            iad=i2e(i2eu)+ijrow+njrow*ijcol
            bc(iho+i)=bc(iho+i)+0.5d0*bc(iad)
            if(ldebug)write(6,*)('j '),bc(iad),0.5d0,iorb,jorb
           end if                                                       6d11s19
          end do
         end do
         jo=jo+nopen
         if(ldebug)write(6,*)bc(iho+i)                                  10d18s20
        end do
        jho=iho-1
        jhc=ihc-1
       end if
       hd=bc(jhc+ibasis(jbasis+2))+bc(jho+ibasis(jbasis+3))
       if(ldebug)write(6,*)('sum of oo and cc parts: '),hd,jhc,jho,
     $      ibasis(jbasis+1),ibasis(jbasis+2),ibasis(jbasis+3)
       hd=hd+shift
       if(ldebug)write(6,*)('with shift '),hd                           7d1s19
       jc=iptr(2,nclop)+nclo*(ibasis(jbasis+2)-1)-1                     6d14s19
       jo=iptr(4,nclop)+nopen*(ibasis(jbasis+3)-1)-1                    6d14s19
       if(ldebug)then                                                   7d1s19
       write(6,*)('mixing between '),(idorb(jc+j),j=1,nclo)
       write(6,*)('and '),(isorb(jo+j),j=1,nopen)
       end if                                                           7d1s19
       do i=1,nclo
        iorb=idorb(jc+i)
        isb=ism(iorb)
        ig=irel(iorb)-1                                                 6d11s19
        njrow=(irefo(isb)*(irefo(isb)+1))/2                             6d11s19
        ijrow=ig+((ig*(ig+1))/2)                                        6d11s19
        do j=1,nopen                                                    6d11s19
         jorb=isorb(jo+j)                                               6d11s19
         jsb=ism(jorb)                                                  6d11s19
         jg=irel(jorb)-1                                                6d11s19
         ijcol=jg+((jg*(jg+1))/2)                                       6d11s19
         i2eu=ifind2(isb,isb,jsb,jsb,icase)                             6d11s19
         iad=i2e(i2eu)+ijrow+njrow*ijcol                                6d11s19
         hd=hd+2d0*bc(iad)                                              6d11s19
         if(ldebug)write(6,*)('eri'),bc(iad),2d0
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
           if(ldebug)write(6,*)('eri'),bc(iad),-1d0
         hd=hd-bc(iad)
        end do                                                          6d11s19
       end do
       if(ldebug)write(6,*)('final spin independent part: '),hd         7d1s19
       hdiag(ihdiag)=hd
       ihdiag=ihdiag+1
      end do
      return
      end
