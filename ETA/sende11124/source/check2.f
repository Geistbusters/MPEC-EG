c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine check2(iff22r,iff22,nff22r,nff22,nfdatr,nfdat,ncsfr,   8d18s22
     $     ncsf,ncsf2r,ncsf2,mff2a,mdon,mdoo,nsymb,norb,ff22r,bc,ibc)   11d14s22
      implicit real*8 (a-h,o-z)                                         8d18s22
      integer*8 iff22r(*),iff22(*),ipack8                               8d18s22
      dimension nff22r(mdoo+1,2,*),nff22(mdoo+1,2,*),nfdatr(5,4,*),     8d18s22
     $     nfdat(5,4,*),ncsfr(*),ncsf(*),ncsf2r(4,*),ncsf2(4,*),        8d18s22
     $     ipack4(2),nl(4),ff22r(*)                                     8d18s22
      equivalence (ipack8,ipack4)                                       8d18s22
      include "common.store"                                            8d18s22
      ier=0                                                             8d18s22
      if(ier.ne.0)then                                                  8d18s22
       write(6,*)('for nff22 ')                                         8d18s22
       write(6,*)('got  :'),(((nff22r(ii,j,isb),ii=mdon+1,mdoo+1),      8d18s22
     $      j=1,2),isb=1,nsymb)                                         8d18s22
       write(6,*)('want  :'),(((nff22(ii,j,isb),ii=mdon+1,mdoo+1),      8d18s22
     $      j=1,2),isb=1,nsymb)                                         8d18s22
       stop 'restart'                                                   8d18s22
      end if                                                            8d18s22
      do ii=1,mdoo+1-mdon                                               8d18s22
       if(ncsfr(ii).ne.ncsf(ii))ier=ier+1                               8d18s22
      end do                                                            8d18s22
      if(ier.ne.0)then                                                  8d18s22
       write(6,*)('for ncsf ')
       write(6,*)('got  :'),(ncsfr(ii),ii=1,mdoo+1-mdon)                8d18s22
       write(6,*)('want :'),(ncsf(ii),ii=1,mdoo+1-mdon)                 8d18s22
       stop 'restart'                                                   8d18s22
      end if                                                            8d18s22
      do ii=1,mdoo+1-mdon                                               8d18s22
       if(ncsf(ii).gt.0)then                                            4d19s23
        do l=1,4                                                         8d18s22
         if(ncsf2r(l,ii).ne.ncsf2(l,ii))ier=ier+1                        8d18s22
        end do                                                           8d18s22
       end if                                                           8d18s22
      end do                                                            8d18s22
      if(ier.ne.0)then                                                  8d18s22
       write(6,*)('for ncsf2 ')                                         8d18s22
       write(6,*)('got  :'),((ncsf2r(l,ii),l=1,4),ii=1,mdoo+1-mdon)     8d18s22
       write(6,*)('want :'),((ncsf2(l,ii),l=1,4),ii=1,mdoo+1-mdon)      8d18s22
       stop 'restart'                                                   8d18s22
      end if                                                            8d18s22
      do isb=1,nsymb                                                    8d18s22
       do nclo=mdon,mdoo                                                8d18s22
        nclop=nclo+1                                                    8d18s22
        if(nff22(nclop,1,isb).gt.0)then                                 8d18s22
         iarg=nclop-mdon                                                8d18s22
         ivcv=nfdat(5,1,isb)+nff22(nclop,2,isb)                         8d18s22
         ivcvr=nfdatr(5,1,isb)+nff22(nclop,2,isb)                       8d18s22
         do if=1,nff22(nclop,1,isb)                                     8d18s22
          ipack8=ibc(ivcv)                                              8d18s22
          itc=ipack4(1)                                                 8d18s22
          ito=ipack4(2)                                                 8d18s22
          ipack8=iff22r(ivcvr)                                          8d18s22
          if(ipack4(1).ne.itc.or.ipack4(2).ne.ito)then                  8d18s22
           write(6,*)('for if,nclo,isb '),if,nclo,isb,(' itc,ito')      8d18s22
           write(6,*)('got  :'),ipack4                                  8d18s22
           write(6,*)('want :'),itc,ito                                 8d18s22
           call dcbit(ipack4,norb,'cgot')                               8d18s22
           call dcbit(ipack4(2),norb,'ogot')                            8d18s22
           call dcbit(itc,norb,'cwant')                                 8d18s22
           call dcbit(ito,norb,'owant')                                 8d18s22
           stop 'restart'                                               8d18s22
          end if                                                        8d18s22
          nspace=ibc(ivcv+1)                                            8d18s22
          if(ibc(ivcv+1).ne.iff22r(ivcvr+1))then                        8d18s22
           write(6,*)('for if,nclo,isb '),if,nclo,isb,(' nspace')       8d18s22
           write(6,*)('got  :'),iff22r(ivcvr+1)                         8d18s22
           write(6,*)('want :'),ibc(ivcv+1)                             8d18s22
           stop 'restart'                                               8d18s22
          end if                                                        8d18s22
          do l=1,4                                                      8d18s22
           nl(l)=ibc(ivcv+1+l)                                          8d18s22
           if(ibc(ivcv+1+l).ne.iff22r(ivcvr+1+l))ier=ier+1              8d18s22
          end do                                                        8d18s22
          if(ier.ne.0)then                                              8d18s22
           write(6,*)('for if,nclo,isb '),if,nclo,isb,(' nl')           8d18s22
           write(6,*)('got  :'),(iff22r(ivcvr+1+l),l=1,4)               8d18s22
           write(6,*)('want :'),(ibc(ivcv+1+l),l=1,4)                   8d18s22
           stop 'restart'                                               8d18s22
          end if                                                        8d18s22
          do l=1,4                                                      8d18s22
           if(nl(l).gt.0)then                                           8d18s22
            iad1=ivcv+ibc(ivcv+5+l)                                      8d18s22
            iad2=iad1+nl(l)                                             8d18s22
            iad1r=ivcvr+iff22r(ivcvr+5+l)                               8d18s22
            iad2r=iad1r+nl(l)                                             8d18s22
            rms=0d0                                                      8d18s22
            do i=0,ncsf2(l,iarg)*nl(l)-1                                 8d18s22
             rms=rms+(bc(iad2+i)-ff22r(iad2r+i))**2                     8d18s22
            end do                                                       8d18s22
            rms=sqrt(rms/dfloat(ncsf2(l,iarg)*nl(l)))                   8d18s22
            if(rms.gt.1d-10)then                                        8d18s22
             write(6,*)('for l,if,nclo,isb '),l,if,nclo,isb,(' vectors')8d18s22
             write(6,*)('got ')                                         8d18s22
             call prntm2(ff22r(iad2r),ncsf2(l,iarg),nl(l),ncsf2(l,iarg))8d18s22
             write(6,*)('want ')                                        8d18s22
             call prntm2(bc(iad2),ncsf2(l,iarg),nl(l),ncsf2(l,iarg))    8d18s22
             stop 'restart'                                               8d18s22
            end if                                                      8d18s22
            do i=0,nl(l)-1                                              8d18s22
             if(ibc(iad1+i).ne.iff22r(iad1r+i))ier=ier+1                8d18s22
            end do                                                      8d18s22
            if(ier.gt.0)then                                            8d18s22
             write(6,*)('for l,if,nclo,isb '),l,if,nclo,isb,
     $           (' pointers')                                          8d18s22
             write(6,*)('got  :'),(iff22r(iad1r+i),i=0,nl(l)-1)         8d18s22
             write(6,*)('want :'),(ibc(iad1+i),i=0,nl(l)-1)             8d18s22
             stop 'restart'                                             8d18s22
            end if                                                      8d18s22
           end if
          end do                                                        8d18s22
          ivcv=ivcv+nspace                                              8d18s22
          ivcvr=ivcvr+nspace                                            8d18s22
         end do                                                         8d18s22
        end if                                                          8d18s22
       end do                                                           8d18s22
      end do                                                            8d18s22
      return                                                            8d18s22
      end                                                               8d18s22
