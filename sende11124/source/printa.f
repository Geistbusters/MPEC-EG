c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine printa(hessa,n1,io1,n2,io2,n3,io3,n4,io4,tmp)
      implicit real*8 (a-h,o-z)
      dimension isym(2,11,4),nsym(4),hessa(*),tmp(*),ps(2)
      data nsym/11,7,4,2/
      data ((isym(k,j,1),k=1,2),j=1,11)/1,1, 2,1, 4,1, 6,1, 9,1,
     $     11,1, 13,1, 14,1, 17,1, 22,1, 23,2/
      data ((isym(k,j,2),k=1,2),j=1,7)/3,1, 7,1, 8,2, 12,1, 18,2,
     $     19,2, 24,1/
      data ((isym(k,j,3),k=1,2),j=1,4)/5,1, 10,1, 16,1, 21,1/
      data ((isym(k,j,4),k=1,2),j=1,2)/15,1, 20,2/
      data ps/1d0,-1d0/
      write(6,*)('printa')
      write(6,*)('starting with ')
      call prntm2(hessa,n1*n2,n3*n4,n1*n2)
      do isa=1,4
       do isb=1,4
        do isc=1,4
         do isd=1,4
          ito=1
          sz=0d0
          nab=0
          ncdx=0
          do ia=1,nsym(isa)
           iam=isym(1,ia,isa)-io4
           if(iam.ge.1.and.iam.le.n4)then
            do ib=1,nsym(isb)
             ibm=isym(1,ib,isb)-io3
             if(ibm.ge.1.and.ibm.le.n3)then
              pab=ps(isym(2,ia,isa))*ps(isym(2,ib,isb))
              ncd=0
              do ic=1,nsym(isc)
               icm=isym(1,ic,isc)-io2
               if(icm.ge.1.and.icm.le.n2)then
                pabc=pab*ps(isym(2,ic,isc))
                do id=1,nsym(isd)
                 idm=isym(1,id,isd)-io1
                 if(idm.ge.1.and.idm.le.n1)then
                  ifrm=idm+n1*(icm-1+n2*(ibm-1+n3*(iam-1)))
                  tmp(ito)=hessa(ifrm)*pabc*ps(isym(2,id,isd))          6d21s22
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
      write(6,*)('returning from printa')
      return
      end
      subroutine printb(hessa,n1,io1,n2,io2,n3,io3,n4,io4,tmp)
      implicit real*8 (a-h,o-z)
      dimension isym(2,11,4),nsym(4),hessa(*),tmp(*),ps(2)
      data nsym/11,7,4,2/
      data ((isym(k,j,1),k=1,2),j=1,11)/1,1, 3,1, 5,1, 7,1, 9,1,
     $     13,1, 14,1, 17,1, 18,2, 22,2, 24,1/
      data ((isym(k,j,2),k=1,2),j=1,7)/2,1, 4,1, 6,2, 12,1, 15,1,
     $     19,2, 23,2/
      data ((isym(k,j,3),k=1,2),j=1,4)/8,1, 10,1, 16,1, 21,1/
      data ((isym(k,j,4),k=1,2),j=1,2)/11,1, 20,2/
      data ps/1d0,-1d0/
      write(6,*)('printb')
      write(6,*)('starting with ')
      call prntm2(hessa,n1*n2,n3*n4,n1*n2)
      do isa=1,4
       do isb=1,4
        do isc=1,4
         do isd=1,4
          ito=1
          sz=0d0
          nab=0
          ncdx=0
          do ia=1,nsym(isa)
           iam=isym(1,ia,isa)-io4
           if(iam.ge.1.and.iam.le.n4)then
            do ib=1,nsym(isb)
             ibm=isym(1,ib,isb)-io3
             if(ibm.ge.1.and.ibm.le.n3)then
              pab=ps(isym(2,ia,isa))*ps(isym(2,ib,isb))
              ncd=0
              do ic=1,nsym(isc)
               icm=isym(1,ic,isc)-io2
               if(icm.ge.1.and.icm.le.n2)then
                pabc=pab*ps(isym(2,ic,isc))
                do id=1,nsym(isd)
                 idm=isym(1,id,isd)-io1
                 if(idm.ge.1.and.idm.le.n1)then
                  ifrm=idm+n1*(icm-1+n2*(ibm-1+n3*(iam-1)))
                  tmp(ito)=hessa(ifrm)*pabc*ps(isym(2,id,isd))          6d21s22
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
      write(6,*)('returning from printb')
      return
      end
      subroutine printba(hessa,n1,io1,n2,io2,n3,io3,n4,io4,tmp)
      implicit real*8 (a-h,o-z)
      dimension isym(2,11,4),nsym(4),hessa(*),tmp(*),ps(2),
     $     jsym(2,11,4)
      data nsym/11,7,4,2/
      data ((jsym(k,j,1),k=1,2),j=1,11)/1,1, 2,1, 5,1, 6,1, 9,1,
     $     11,1, 13,1, 14,1, 18,1, 22,1, 23,1/
      data ((jsym(k,j,2),k=1,2),j=1,7)/3,1, 7,2, 8,2, 12,1, 17,2,
     $     19,2, 24,1/
      data ((jsym(k,j,3),k=1,2),j=1,4)/4,1, 10,1, 16,1, 21,1/
      data ((jsym(k,j,4),k=1,2),j=1,2)/15,1, 20,2/
      data ((isym(k,j,1),k=1,2),j=1,11)/1,1, 3,1, 5,1, 7,1, 9,1,
     $     13,1, 14,1, 17,1, 18,2, 22,2, 24,1/
      data ((isym(k,j,2),k=1,2),j=1,7)/2,1, 4,1, 6,2, 12,1, 15,1,
     $     19,2, 23,2/
      data ((isym(k,j,3),k=1,2),j=1,4)/8,1, 10,1, 16,1, 21,1/
      data ((isym(k,j,4),k=1,2),j=1,2)/11,1, 20,2/
      data ps/1d0,-1d0/
      write(6,*)('printba')
      write(6,*)('starting with ')
      call prntm2(hessa,n1*n2,n3*n4,n1*n2)
      do isa=1,4
       do isb=1,4
        do isc=1,4
         do isd=1,4
          ito=1
          sz=0d0
          nab=0
          ncdx=0
          do ia=1,nsym(isa)
           iam=jsym(1,ia,isa)-io4
           if(iam.ge.1.and.iam.le.n4)then
            do ib=1,nsym(isb)
             ibm=jsym(1,ib,isb)-io3
             if(ibm.ge.1.and.ibm.le.n3)then
              pab=ps(jsym(2,ia,isa))*ps(jsym(2,ib,isb))
              ncd=0
              do ic=1,nsym(isc)
               icm=isym(1,ic,isc)-io2
               if(icm.ge.1.and.icm.le.n2)then
                pabc=pab*ps(isym(2,ic,isc))
                do id=1,nsym(isd)
                 idm=isym(1,id,isd)-io1
                 if(idm.ge.1.and.idm.le.n1)then
                  ifrm=idm+n1*(icm-1+n2*(ibm-1+n3*(iam-1)))
                  tmp(ito)=hessa(ifrm)*pabc*ps(isym(2,id,isd))          6d21s22
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
            write(6,*)('times2! ')
            do i=1,ncd*nab
             tmp(i)=tmp(i)*2d0
            end do
            call prntm2(tmp,ncd,nab,ncd)
           end if
          end if
         end do
        end do
       end do
      end do
      write(6,*)('returning from printba')
      return
      end
      subroutine printc(hessa,tmp)
      implicit real*8 (a-h,o-z)
      dimension isym(2,11,4),nsym(4),hessa(24,24),tmp(*),ps(2),
     $     jsym(2,11,4)
      data nsym/11,7,4,2/
      data ((isym(k,j,1),k=1,2),j=1,11)/1,1, 2,1, 3,1, 6,1, 9,1,
     $     10,1, 11,1, 15,1, 16,1, 21,1, 22,2/
      data ((isym(k,j,2),k=1,2),j=1,7)/4,1, 7,1, 13,1, 17,2, 18,2,
     $     19,1, 24,2/
      data ((isym(k,j,3),k=1,2),j=1,4)/5,1, 8,1, 14,1, 20,1/
      data ((isym(k,j,4),k=1,2),j=1,2)/12,1, 23,2/
      data ps/1d0,-1d0/
      write(6,*)('printc')
      write(6,*)('starting with ')
      call prntm2(hessa,24,24,24)
      srh=sqrt(0.5d0)
      write(6,*)('take rows from cao to sao')
      do i=1,14
       write(6,*)('copy row '),i
       do j=1,24
        ij=i+24*(j-1)
        tmp(ij)=hessa(i,j)
       end do
      end do
      do i=15,16
       ip=i+2
       write(6,*)('combine rows '),i,ip
       do j=1,24
        ij=i+24*(j-1)
        ipj=ip+24*(j-1)
        tmp(ij)=srh*(hessa(i,j)+hessa(ip,j))
        tmp(ipj)=srh*(hessa(i,j)-hessa(ip,j))
       end do
      end do
      do i=19,21
       ip=i+3
       write(6,*)('combine rows '),i,ip
       do j=1,24
        ij=i+24*(j-1)
        ipj=ip+24*(j-1)
        tmp(ij)=srh*(hessa(i,j)+hessa(ip,j))
        tmp(ipj)=srh*(hessa(i,j)-hessa(ip,j))
       end do
      end do
      write(6,*)('new matrix ')
      call prntm2(tmp,24,24,24)
      ioff=24*24
      write(6,*)('take cols from cao to sao')
      do i=1,14
       write(6,*)('copy col '),i
       do j=1,24
        ji=j+24*(i-1)
        tmp(ji+ioff)=tmp(ji)
       end do
      end do
      do i=15,16
       ip=i+2
       write(6,*)('combine cols '),i,ip
       do j=1,24
        ji=j+24*(i-1)
        jip=j+24*(ip-1)
        tmp(ji+ioff)=srh*(tmp(ji)+tmp(jip))
        tmp(jip+ioff)=srh*(tmp(ji)-tmp(jip))
       end do
      end do
      do i=19,21
       ip=i+3
       write(6,*)('combine cols '),i,ip
       do j=1,24
        ji=j+24*(i-1)
        jip=j+24*(ip-1)
        tmp(ji+ioff)=srh*(tmp(ji)+tmp(jip))
        tmp(jip+ioff)=srh*(tmp(ji)-tmp(jip))
       end do
      end do
      write(6,*)('new matrix ')
      call prntm2(tmp(ioff+1),24,24,24)
      do isa=1,4
       do isc=1,4
        ito=0
        sz=0d0
        do ia=1,nsym(isa)
         iam=isym(1,ia,isa)
         pab=ps(isym(2,ia,isa))
         do ic=1,nsym(isc)
          icm=isym(1,ic,isc)
          pac=pab*ps(isym(2,ic,isc))*2d0
          ifrm=ioff+icm+24*(iam-1)                                           2d27s23
          ito=ito+1
          tmp(ito)=tmp(ifrm)*pac                                        2d27s23
          sz=sz+tmp(ito)**2
         end do
        end do
        if(ito.gt.0)then
         sz=sqrt(sz/dfloat(ito))
         if(sz.gt.1d-5)then
          write(6,*)('times 2')
          write(6,*)('symmetry type '),isc,isa
          call prntm2(tmp,nsym(isc),nsym(isa),nsym(isc))
         end if
        end if
       end do
      end do
      write(6,*)('returning from printc')
      return
      end
