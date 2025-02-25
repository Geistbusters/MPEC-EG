c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine lzziicsfcas(hdiag,nfcn,mdon,mdoo,ibasis,               8d2s22
     $     ixlzz,shift,nec,islz,multh,ism,irel,irefo,idoubo,npass,      8d2s22
     $     iptrbit,mysym,norb,bc,ibc)                                   11d10s22
c
c     compute spin independent part of lzlzii                           12d24s19
c
      implicit real*8 (a-h,o-z)
      logical ldebug                                                    7d1s19
      dimension hdiag(nfcn),ixlzz(8,*),ibasis(*),multh(8,8),            8d2s22
     $     ism(*),irel(*),irefo(*),idoubo(*),islz(*),iptrbit(2,mdoo+1,*)8d2s22
      include "common.store"
      ldebug=.false.                                                    7d15s19
      ione=npass+1                                                      12d31s19
      if(ldebug)write(6,*)('hi, i''m lzziicsf'),nfcn,(islz(i),i=1,npass)12d31s19
     $     ,loc(ixlzz)                                                  12d31s19
      ihc=ibcoff                                                        8d2s22
      if(ldebug)then
       do isb=1,4
        if(irefo(isb).gt.0)then                                         12d31s19
         write(6,*)('ixlzz addresses: '),isb,(ixlzz(isb,j),j=1,npass+1)  12d31s19
         write(6,*)('one electron part: ')
         call prntm2(bc(ixlzz(isb,ione)),irefo(isb),irefo(isb),
     $       irefo(isb))
         do ipass=1,npass                                                12d31s19
          jsb=multh(islz(ipass),isb)
          if(irefo(jsb).gt.0)then                                       12d31s19
           write(6,*)('for two electron matrix '),ipass,jsb
           call prntm2(bc(ixlzz(isb,ipass)),irefo(isb),irefo(jsb),
     $          irefo(isb))
          end if
         end do                                                          12d31s19
        end if
       end do
      end if                                                            7d15s19
      nclo=-1
      ihdiag=1
      do if=1,nfcn
       jbasis=3*(if-1)
       mclo=ibasis(jbasis+1)                                            8d2s22
       mclop=mclo+1                                                     8d2s22
       nopen=nec-2*mclo                                                 8d2s22
       nclo=mclo
       nclop=mclop
       jjc=iptrbit(1,mclop,mysym)+ibasis(jbasis+2)-1                    8d2s22
       jjo=iptrbit(2,mclop,mysym)+ibasis(jbasis+3)-1                    8d2s22
       xc=0d0                                                           8d2s22
       xo=0d0                                                           8d2s22
       do i=1,norb                                                      8d2s22
        isb=ism(i)                                                      8d2s22
        ig=irel(i)-1                                                    8d2s22
        ndm=irefo(isb)                                                  8d2s22
        ih=ixlzz(isb,ione)+ig*(ndm+1)                                   8d2s22
        if(btest(ibc(jjc),i))then                                       8d2s22
         orig=xc
         xc=xc+2d0*bc(ih)                                               8d2s22
         do j=1,norb                                                    8d2s22
          if(btest(ibc(jjc),j))then                                     8d2s22
           jsb=ism(j)                                                   8d2s22
           jg=irel(j)-1                                                 8d2s22
           do ipass=1,npass                                             12d31s19
            if(multh(isb,jsb).eq.islz(ipass))then                       12d31s19
             ndm=irefo(isb)                                             12d31s19
             iad=ixlzz(isb,ipass)+ig+ndm*jg                             12d31s19
             orig=xc
             xc=xc-2d0*bc(iad)*bc(iad)                                  8d2s22
            end if                                                       12d24s19
           end do                                                       12d31s19
          end if                                                        8d2s22
         end do                                                         8d2s22
        else if(btest(ibc(jjo),i))then                                  8d2s22
         orig=xo
         xo=xo+bc(ih)                                                   8d2s22
        end if                                                          8d2s22
       end do                                                           8d2s22
       xd=xc+xo                                                         8d2s22
       hd=xd                                                            8d2s22
       if(ldebug)write(6,*)('sum of oo and cc parts: '),hd              7d1s19
       jjc=iptrbit(1,nclop,mysym)+ibasis(jbasis+2)-1                    8d2s22
       jjo=iptrbit(2,nclop,mysym)+ibasis(jbasis+3)-1                    8d2s22
       do i=1,norb                                                      8d2s22
        if(btest(ibc(jjc),i))then                                       8d2s22
         isb=ism(i)                                                     8d2s22
         ig=irel(i)-1                                                   8d2s22
         do j=1,norb                                                     8d2s22
          if(btest(ibc(jjo),j))then                                      8d2s22
           jsb=ism(j)                                                    8d2s22
           jg=irel(j)-1                                                  8d2s22
           do ipass=1,npass                                               12d31s19
            if(multh(isb,jsb).eq.islz(ipass))then                         12d31s19
             ndm=irefo(isb)                                               12d31s19
             iad=ixlzz(isb,ipass)+ig+ndm*jg                               12d31s19
             hd=hd-2d0*bc(iad)*bc(iad)                                     12d24s19
            end if                                                         12d24s19
           end do                                                         12d31s19
          end if                                                         8d2s22
         end do                                                          6d11s19
        end if                                                          8d2s22
       end do
       if(ldebug)write(6,*)('final spin independent part: '),hd         7d1s19
       hdiag(ihdiag)=hd
       ihdiag=ihdiag+1
      end do
      if(ldebug)then
       call dws_synca
       call dws_finalize
       stop
      end if
      return
      end
