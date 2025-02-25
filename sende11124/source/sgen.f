c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine sgen(ia,id2,id1,nox,nex,n1,iorb)
      implicit real*8 (a-h,o-z)
c
c     generate fci orbital occupancies
c     no is number of orbitals
c     ne is number of electrons
c
      integer*1 ia(id1,id2),iorb(nex,*)
      iorb(1,1)=nox
      no=nox
      ne=nex
c
c     add electrons one at a time into orbitals after the previous electron
c
c     first electron
c
      ntop=no-ne+1
      do i=1,ntop
       do j=1,no
        ia(j,i)=0
       end do
       ia(i,i)=1
      end do
      n1=ntop
c
c     remaining electrons
c
      do ie=2,ne
       ntopn=ntop+1
       in=n1+1
       do ix=1,n1
        do i=ntop,1,-1
         if(ia(i,ix).eq.1)then
          do j=i+2,ntopn
           if(in.gt.id2)stop 'id2'
           do k=1,no
            ia(k,in)=ia(k,ix)
           end do
           ia(j,in)=1
           in=in+1
          end do
          ia(i+1,ix)=1
          go to 3
         end if
        end do
    3   continue
       end do
       n1=in-1
       ntop=ntopn
      end do
      do i=1,n1
       k=1
       do j=1,nox
        if(ia(j,i).eq.1)then
         if(k.gt.nex)stop 'sgen'
         iorb(k,i)=j
         k=k+1
        end if
       end do
      end do
      return
      end
