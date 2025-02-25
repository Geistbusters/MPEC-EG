c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine deramat(i4od,ionexd,iamat,noc,nvirtc,ipuse,ih0der,
     $     multh,idwsdeb,isblkder,isblkxder,nsblkder,nsblkxder,bc,ibc)  11d14s22
      implicit real*8 (a-h,o-z)
c
c     given derivative integrals, compute der of amat.
c     for reference, amat is half the gradient.
c
      include "common.hf"
      dimension i4od(1),ionexd(1),iamat(1),noc(1),nvirtc(1),multh(8,8),
     $     isblkder(4,idbk),isblkxder(4,idbk)                           8d3s16
      include "common.store"
      if(ipuse.eq.1)then
       jh0der=ih0der
       do isb=1,nsymb
        iamat(isb)=ibcoff
        ibcoff=iamat(isb)+nvirtc(isb)*noc(isb)
        igoal=iamat(1)
        nh=noc(isb)+nvirtc(isb)
        if(idwsdeb.gt.1)then
        end if
        call enough('deramat.  1',bc,ibc)
        nh=nvirtc(isb)+noc(isb)
        do icol=0,noc(isb)-1
         do irow=0,nvirtc(isb)-1
          iadh=jh0der+noc(isb)+irow+nh*icol                             1d3s17
          iada=iamat(isb)+irow+nvirtc(isb)*icol
          bc(iada)=bc(iadh)*2d0
         end do
        end do
        if(idwsdeb.gt.1)then
         write(6,*)('der of amat after dH: '),isb,iamat(isb)
         call prntm2(bc(iamat(isb)),nvirtc(isb),noc(isb),nvirtc(isb))
         if(nsymb.eq.1)call printa(bc(iamat(isb)),nvirtc(isb),noc(isb),
     $        1,0,noc(isb),0,1,0,bc(ibcoff))
        end if
        do is=1,nsblkxder
         if(isblkxder(1,is).eq.isblkxder(2,is).and.
     $        isblkxder(3,is).eq.isb)then
          n2=noc(isblkxder(1,is))*noc(isblkxder(1,is))
          np=noc(isblkxder(1,is))+1
          do ir=0,nvirtc(isb)-1
           do m=0,noc(isb)-1
            iadx=ionexd(is)+n2*(m+noc(isb)*ir)
            sum=0d0
            do il=0,noc(isblkxder(1,is))-1
             iad1=iadx+il*np
             sum=sum+bc(iad1)
            end do
            iada=iamat(isb)+ir+nvirtc(isb)*m
            bc(iada)=bc(iada)+4d0*sum
           end do
          end do
         end if
         if(isblkxder(1,is).eq.isb.and.
     $          isblkxder(2,is).eq.isblkxder(3,is))then
          n3=noc(isb)*noc(isblkxder(2,is))*noc(isblkxder(2,is))
          np=(noc(isblkxder(2,is))+1)*noc(isb)
          do ir=0,nvirtc(isb)-1
           do m=0,noc(isb)-1
            iadx=ionexd(is)+m+n3*ir
            sum=0d0
            do il=0,noc(isblkxder(2,is))-1
             iad1=iadx+np*il
             sum=sum+bc(iad1)
            end do
            iada=iamat(isb)+ir+nvirtc(isb)*m
            if(isblkxder(1,is).ne.isblkxder(2,is))sum=sum*2d0
            bc(iada)=bc(iada)-sum
           end do
          end do
         end if
         if(isblkxder(2,is).eq.isb.and.
     $          isblkxder(1,is).eq.isblkxder(3,is))then
          n2=noc(isb)*noc(isblkxder(1,is))
          np=n2+1
          do ir=0,nvirtc(isb)-1
           do m=0,noc(isb)-1
            iadx=ionexd(is)+noc(isblkxder(1,is))*(m+n2*ir)
            sum=0d0
            do il=0,noc(isblkxder(1,is))-1
             iad1=iadx+np*il
             sum=sum+bc(iad1)
            end do
            iada=iamat(isb)+ir+nvirtc(isb)*m
            if(isblkxder(1,is).ne.isblkxder(2,is))sum=sum*2d0
            bc(iada)=bc(iada)-sum
           end do
          end do
         end if
        end do
        if(idwsdeb.gt.10)then
         write(6,*)('der of amat for symmetry block '),isb,iamat(isb)
         call prntm2(bc(iamat(isb)),nvirtc(isb),noc(isb),nvirtc(isb))   1d3s17
         if(nsymb.eq.1)call printa(bc(iamat(1)),nvirtc,noc,1,0,noc,0,1,
     $        0,bc(ibcoff))
        end if
        jh0der=jh0der+nh*nh
       end do
      else
       jh0der=ih0der
       do isb=1,nsymb                                                   5d2s22
        isk=multh(ipuse,isb)
        iamat(isb)=ibcoff                                               5d2s22
        ibcoff=iamat(isb)+nvirtc(isb)*noc(isk)                          5d2s22
       end do                                                           5d2s22
       call enough('deramat.  2',bc,ibc)
       do isb=1,nsymb
        isk=multh(ipuse,isb)
        if(isk.gt.isb)then
         nk=nvirtc(isk)+noc(isk)
         nb=nvirtc(isb)+noc(isb)
         if(idwsdeb.gt.1)then
          write(6,*)('der of h0 ')
          call prntm2(bc(jh0der),nb,nk,nb)
         end if
         call enough('deramat.  3',bc,ibc)
         do icol=0,noc(isk)-1
          do irow=0,nvirtc(isb)-1
           iadh=jh0der+irow+noc(isb)+nb*icol
           iada=iamat(isb)+irow+nvirtc(isb)*icol
           bc(iada)=2d0*bc(iadh)
          end do
         end do
         do icol=0,nvirtc(isk)-1
          do irow=0,noc(isb)-1
           iadh=jh0der+irow+nb*(icol+noc(isk))
           iada=iamat(isk)+icol+nvirtc(isk)*irow
           bc(iada)=2d0*bc(iadh)
          end do
         end do
         jh0der=jh0der+nb*nk
        end if
        if(idwsdeb.gt.1)then
         write(6,*)('der of amat after h part '),isb,isk
         write(6,*)('bra first')
         call prntm2(bc(iamat(isb)),nvirtc(isb),noc(isk),nvirtc(isb))
         write(6,*)('ket first')
         call prntm2(bc(iamat(isk)),nvirtc(isk),noc(isb),nvirtc(isk))
        end if
       end do
       do isb=1,nsymb
        isk=multh(ipuse,isb)
        do is=1,nsblkxder
         if(isblkxder(1,is).eq.isblkxder(2,is).and.
     $        isblkxder(3,is).eq.isk)then
          n2=noc(isblkxder(1,is))**2
          np=1+noc(isblkxder(1,is))
          do m=0,noc(isk)-1
           do ir=0,nvirtc(isb)-1
            sum=0d0
            iadx=ionexd(is)+n2*(m+noc(isk)*ir)
            do il=0,noc(isblkxder(1,is))-1
             iadu=iadx+il*np
             sum=sum+bc(iadu)
            end do
            iada=iamat(isb)+ir+nvirtc(isb)*m
            bc(iada)=bc(iada)+4d0*sum
           end do
          end do
         end if
         if(isblkxder(2,is).eq.isblkxder(3,is).and.
     $        isblkxder(1,is).eq.isk)then
          n3=noc(isk)*noc(isblkxder(2,is))*noc(isblkxder(2,is))
          np=noc(isk)*(1+noc(isblkxder(2,is)))
          do m=0,noc(isk)-1
           do ir=0,nvirtc(isb)-1
            sum=0d0
            iadx=ionexd(is)+m+n3*ir
            do il=0,noc(isblkxder(2,is))-1
             iadu=iadx+il*np
             sum=sum+bc(iadu)
            end do
            if(isblkxder(1,is).ne.isblkxder(2,is))sum=sum*2d0
            iada=iamat(isb)+ir+nvirtc(isb)*m
            bc(iada)=bc(iada)-sum
           end do
          end do
         end if
         if(isblkxder(1,is).eq.isblkxder(3,is).and.
     $        isblkxder(2,is).eq.isk)then
          n2=noc(isk)*noc(isblkxder(1,is))
          np=1+noc(isk)*noc(isblkxder(1,is))
          do m=0,noc(isk)-1
           do ir=0,nvirtc(isb)-1
            sum=0d0
            iadx=ionexd(is)+noc(isblkxder(1,is))*(m+n2*ir)
            do il=0,noc(isblkxder(1,is))-1
             iadu=iadx+il*np
             sum=sum+bc(iadu)
            end do
            if(isblkxder(1,is).ne.isblkxder(2,is))sum=sum*2d0
            iada=iamat(isb)+ir+nvirtc(isb)*m
            bc(iada)=bc(iada)-sum
           end do
          end do
         end if
        end do
        if(idwsdeb.gt.10)then
        write(6,*)('for symmetry block '),isb,isk
        write(6,*)('der of amat: ')
        write(6,*)('bra first')
        call prntm2(bc(iamat(isb)),nb,noc(isk),nb)
        write(6,*)('ket first')
        call prntm2(bc(iamat(isk)),nk,noc(isb),nk)
        end if
       end do
      end if
      return
      end
