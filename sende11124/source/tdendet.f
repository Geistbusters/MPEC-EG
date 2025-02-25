c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine tdendet(vec0,n0,is0,vecx,nx,isx,nadet,nbdet,iacto,     3d17s23
     $     jden1i,iaorb,nalpha,iborb,nbeta,multh,nsymb,norb,ism,irel,   3d17s23
     $     bc,ibc,igoal,ioffdeta,ioffdetb)                                                      3d17
      implicit real*8 (a-h,o-z)
      integer*1 iaorb,iborb                                             3d17s23
      integer*8 gandc,ism,irel                                          3d17s23
      dimension vec0(n0,*),is0(*),vecx(nx,*),isx(*),nadet(*),nbdet(*),  3d17s23
     $     iacto(*),jden1i(*),iaorb(nalpha,*),iborb(nbeta,*),multh(8,8),3d17s23
     $     iaoff(8),iboff(8),ism(*),irel(*),ivoff(8,2),ioffdeta(*),     3d20s23
     $     ioffdetb(*)                                                  3d20s23
      include "common.store"
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      i0=0                                                              3d17s23
      ix=0                                                              3d17s23
      do isb=1,nsymb                                                    3d17s23
       ivoff(isb,1)=i0                                                  3d17s23
       ivoff(isb,2)=ix                                                  3d17s23
       call ilimts(nadet(isb),1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d17s23
     $      i2e)                                                        3d17s23
       jsb0=multh(isb,is0(1))                                           3d21s23
       jsbx=multh(isb,isx(1))                                           3d21s23
       if(ih.ge.il)then                                                 3d17s23
        ibb0=ibcoff
        ibbx=ibb0+nbdet(jsb0)
        ibcoff=ibbx+nbdet(jsbx)
        call enough('testdet.bb',bc,ibc)
        jbb0=ibb0-1                                                      3d17s23
        jbbx=ibbx-1                                                      3d17s23
        do ib=1,nbdet(jsb0)                                              3d17s23
         ibc(jbb0+ib)=0                                                  3d17s23
         do i=1,nbeta                                                    3d17s23
          ibc(jbb0+ib)=ibset(ibc(jbb0+ib),iborb(i,ioffdetb(jsb0)+ib))   3d20s23
         end do                                                          3d17s23
        end do                                                           3d17s23
        do ib=1,nbdet(jsbx)                                              3d17s23
         ibc(jbbx+ib)=0                                                  3d17s23
         do i=1,nbeta                                                    3d17s23
          ibc(jbbx+ib)=ibset(ibc(jbbx+ib),iborb(i,ioffdetb(jsbx)+ib))   3d20s23
         end do                                                          3d17s23
        end do                                                           3d17s23
        do ib=1,nbdet(jsb0)                                              3d17s23
         do jb=1,nbdet(jsbx)                                             3d17s23
          gandc=ieor(ibc(jbbx+jb),ibc(jbb0+ib))                          3d17s23
          ndif=popcnt(gandc)                                             3d17s23
          if(ndif.eq.2)then                                              3d17s23
           ii=0                                                          3d17s23
           do i=1,norb                                                   3d17s23
            if(btest(gandc,i))then                                       3d17s23
             if(btest(ibc(jbbx+jb),i))then                               3d17s23
              iorbx=i                                                     3d17s23
              do j=i+1,norb                                              3d17s23
               if(btest(ibc(jbbx+jb),j))ii=ii+1                          3d17s23
              end do                                                     3d17s23
             else                                                        3d17s23
              iorb0=i                                                     3d17s23
              do j=i+1,norb                                              3d17s23
               if(btest(ibc(jbb0+ib),j))ii=ii+1                          3d17s23
              end do                                                     3d17s23
             end if                                                      3d17s23
            end if                                                       3d17s23
           end do                                                        3d17s23
           if(mod(ii,2).eq.0)then                                        3d17s23
            phs=1d0                                                      3d17s23
           else                                                          3d17s23
            phs=-1d0                                                     3d17s23
           end if                                                        3d17s23
           isb0=ism(iorb0)                                                 3d17s23
           isbx=ism(iorbx)                                                 3d17s23
           iad=jden1i(isb0)+irel(iorb0)-1+iacto(isb0)*(irel(iorbx)-1)    3d17s23
           nad=iacto(isb0)*iacto(isbx)                                   3d17s23
           do irx=1,isx(4)                                               3d17s23
            do ir0=1,is0(4)
             sum=0d0                                                    3d17s23
             do ia=il,ih                                                3d17s23
              iad0=i0+ib+nbdet(jsb0)*(ia-1)                             3d17s23
              iadx=ix+jb+nbdet(jsbx)*(ia-1)                             3d17s23
              sum=sum+vec0(iad0,ir0)*vecx(iadx,irx)                     3d17s23
             end do                                                     3d17s23
             bc(iad)=bc(iad)+sum*phs                                    3d17s23
             iad=iad+nad                                                3d17s23
            end do                                                      3d17s23
           end do                                                       3d17s23
          end if                                                         3d17s23
         end do                                                          3d17s23
        end do                                                           3d17s23
        ibcoff=ibb0                                                     3d17s23
       end if                                                           3d17s23
       i0=i0+nbdet(jsb0)*nadet(isb)                                     3d17s23
       ix=ix+nbdet(jsbx)*nadet(isb)                                     3d17s23
      end do                                                            3d17s23
      do jsb=1,nsymb                                                    3d17s23
       call ilimts(nbdet(jsb),1,mynprocg,mynowprog,il,ih,i1s,i1e,i2s,   3d17s23
     $      i2e)                                                        3d17s23
       isb0=multh(jsb,is0(1))                                           3d21s23
       isbx=multh(jsb,isx(1))                                           3d21s23
       if(ih.ge.il)then                                                 3d17s23
        iaa0=ibcoff                                                     3d17s23
        iaax=iaa0+nadet(isb0)                                           3d17s23
        ibcoff=iaax+nadet(isbx)                                         3d17s23
        call enough('testdet.aa',bc,ibc)                                3d17s23
        jaa0=iaa0-1                                                     3d17s23
        jaax=iaax-1                                                     3d17s23
        do ia=1,nadet(isb0)                                             3d17s23
         ibc(jaa0+ia)=0                                                 3d17s23
         do i=1,nalpha                                                  3d17s23
          ibc(jaa0+ia)=ibset(ibc(jaa0+ia),iaorb(i,ioffdeta(isb0)+ia))   3d20s23
         end do                                                         3d17s23
        end do                                                          3d17s23
        do ia=1,nadet(isbx)                                             3d17s23
         ibc(jaax+ia)=0                                                 3d17s23
         do i=1,nalpha                                                  3d17s23
          ibc(jaax+ia)=ibset(ibc(jaax+ia),iaorb(i,ioffdeta(isbx)+ia))   3d20s23
         end do                                                         3d17s23
        end do                                                          3d17s23
        do ia=1,nadet(isb0)                                             3d17s23
         do ja=1,nadet(isbx)                                            3d17s23
          gandc=ieor(ibc(jaax+ja),ibc(jaa0+ia))                         3d17s23
          ndif=popcnt(gandc)                                            3d17s23
          if(ndif.eq.2)then                                             3d17s23
           ii=0                                                         3d17s23
           do i=1,norb                                                  3d17s23
            if(btest(gandc,i))then                                      3d17s23
             if(btest(ibc(jaax+ja),i))then                              3d17s23
              jorbx=i                                                   3d17s23
              do j=i+1,norb                                             3d17s23
               if(btest(ibc(jaax+ja),j))ii=ii+1                         3d17s23
              end do                                                    3d17s23
             else                                                       3d17s23
              jorb0=i                                                   3d17s23
              do j=i+1,norb                                             3d17s23
               if(btest(ibc(jaa0+ia),j))ii=ii+1                         3d17s23
              end do                                                    3d17s23
             end if                                                     3d17s23
            end if                                                      3d17s23
           end do                                                       3d17s23
           if(mod(ii,2).eq.0)then                                       3d17s23
            phs=1d0                                                     3d17s23
           else                                                         3d17s23
            phs=-1d0                                                    3d17s23
           end if                                                       3d17s23
           jsb0=ism(jorb0)                                                 3d17s23
           jsbx=ism(jorbx)                                                 3d17s23
           iad=jden1i(jsb0)+irel(jorb0)-1+iacto(jsb0)*(irel(jorbx)-1)    3d17s23
           nad=iacto(jsb0)*iacto(jsbx)                                   3d17s23
           do irx=1,isx(4)                                               3d17s23
            do ir0=1,is0(4)
             sum=0d0                                                    3d17s23
             do ib=il,ih                                                3d17s23
              iad0=ivoff(isb0,1)+ib+nbdet(jsb)*(ia-1)                             3d17s23
              iadx=ivoff(isbx,2)+ib+nbdet(jsb)*(ja-1)                             3d17s23
              sum=sum+vec0(iad0,ir0)*vecx(iadx,irx)                     3d17s23
             end do                                                     3d17s23
             bc(iad)=bc(iad)+sum*phs                                    3d17s23
             iad=iad+nad                                                3d17s23
            end do                                                      3d17s23
           end do                                                       3d17s23
          end if                                                        3d17s23
         end do                                                         3d17s23
        end do                                                          3d17s23
        ibcoff=iaa0
       end if                                                           3d17s23
      end do                                                            3d17s23
      do isb=1,nsymb                                                    3d17s23
       jsb=multh(isb,multh(is0(1),isx(1)))                              3d17s23
       jden1i(isb)=jden1i(isb)+is0(4)*isx(4)*iacto(isb)*iacto(jsb)      3d17s23
      end do                                                            3d17s23
      return                                                            3d17s23
      end                                                               3d17s23
      subroutine examinev(vec,nvec,nroot,nadet,nbdet,multh,mysym,
     $     iaorb,nalpha,iborb,nbeta,nsymb,ioffdeta,ioffdetb,ism,irel)            3d20s23
      implicit real*8 (a-h,o-z)
      integer*1 iaorb(nalpha,*),iborb(nbeta,*)
      integer*8 ism(*),irel(*)
      dimension vec(nvec,*),nadet(*),nbdet(*),multh(8,8),ioffdeta(*),   3d20s23
     $     ioffdetb(*)
      write(6,*)('ioffdeta: '),(ioffdeta(isb),isb=1,nsymb)
      write(6,*)('ioffdetb: '),(ioffdetb(isb),isb=1,nsymb)
      write(6,*)('irel: '),(irel(j),j=1,6)
      write(6,*)('ism : '),(ism(j),j=1,6)
      ivoff=0
      do isb=1,nsymb
       jsb=multh(isb,mysym)
       write(6,*)nadet(isb),jsb,nbdet(jsb)
       iaoff=ioffdeta(isb)
       iboff=ioffdetb(jsb)
       do ia=1,nadet(isb)
        isga=1
        do j=1,nalpha
         isga=multh(isga,ism(iaorb(j,iaoff+ia)))
        end do
        do ib=1,nbdet(jsb)
         isgb=1
         do j=1,nbeta
          isgb=multh(isgb,ism(iborb(j,iboff+ib)))
         end do
         write(6,*)ivoff+ib,(vec(ivoff+ib,ir),ir=1,nroot),
     $        (iaorb(j,iaoff+ia),j=1,nalpha),(':'),
     $        (iborb(j,iboff+ib),j=1,nbeta),('swa: '),isb,isga,
     $        ('swb: '),jsb,isgb
        end do
        ivoff=ivoff+nbdet(jsb)
       end do
      end do
      return
      end
