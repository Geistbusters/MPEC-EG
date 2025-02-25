c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine printat(hessa,n1,tmp)
      implicit real*8 (a-h,o-z)
      dimension isym(2,11,7),nsym(7),hessa(*),tmp(*),ps(2)
      data nsym/8,3,3,1,5,2,2/
      data ((isym(k,j,1),k=1,2),j=1,8)/1,1, 2,1, 6,1, 9,1, 12,1, 16,1,
     $     20,1, 24,1/
      data ((isym(k,j,2),k=1,2),j=1,3)/4,1, 10,1, 18,1/
      data ((isym(k,j,3),k=1,2),j=1,3)/5,1, 11,1, 19,1/
      data ((isym(k,j,4),k=1,2),j=1,1)/17,1/
      data ((isym(k,j,5),k=1,2),j=1,5)/3,1, 7,1, 8,1, 13,1, 23,1/
      data ((isym(k,j,6),k=1,2),j=1,2)/14,1, 21,1/
      data ((isym(k,j,7),k=1,2),j=1,2)/15,1, 22,1/
      data ps/1d0,-1d0/
      write(6,*)('starting with ')
      nn=(n1*(n1+1))/2
      call prntm2(hessa,nn,nn,nn)
      do isa=1,7
       do isb=1,7
        do isc=1,7
         do isd=1,7
          ito=1
          sz=0d0
          nab=0
          ncdx=0
          do ia=1,nsym(isa)
           iam=isym(1,ia,isa)-1
           if(iam.ge.1.and.iam.le.n1)then
            do ib=1,nsym(isb)
             ibm=isym(1,ib,isb)-1
             if(ibm.ge.1.and.ibm.le.n1)then
              pab=ps(isym(2,ia,isa))*ps(isym(2,ib,isb))
              ncd=0
              do ic=1,nsym(isc)
               icm=isym(1,ic,isc)-1
               if(icm.ge.1.and.icm.le.n1)then
                pabc=pab*ps(isym(2,ic,isc))
                do id=1,nsym(isd)
                 idm=isym(1,id,isd)-1
                 if(idm.ge.1.and.idm.le.n1)then
                  ifrm=idm+n1*(icm-1+n2*(ibm-1+n3*(iam-1)))
                  ix=max(idm,icm)
                  in=min(idm,icm)
                  irow=((ix*(ix-1))/2)+in
                  ix=max(iam,ibm)
                  in=min(iam,ibm)
                  icol=((ix*(ix-1))/2)+in-1
                  ifrm=irow+nn*icol
                  tmp(ito)=hessa(ifrm)*pabc*ps(isym(2,id,isd))
                  sz=sz+tmp(ito)**2
                  ito=ito+1
                  ncd=ncd+1
                  ncdx=max(ncdx,ncd)
                 end if
                end do
               end if
              end do
              nab=nab+1
             end if
            end do
           end if
          end do
          if(ito.gt.0)then
           sz=sqrt(sz/dfloat(ito))
           if(sz.gt.1d-5)then
            write(6,*)('symmetry type '),isd,isc,isb,isa,ifrm,
     $           hessa(ifrm)
            call prntm2(tmp,ncd,nab,ncd)
           end if
          end if
         end do
        end do
       end do
      end do
      return
      end
