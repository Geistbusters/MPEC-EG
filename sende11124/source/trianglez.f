c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine trianglez(i4o,ionex,noc,nvirt,nsblk,isblk,             5d9s22
     $     nsblkx,isblkx,idwsdeb,dob,bc,ibc)                            11d10s22
      implicit real*8 (a-h,o-z)                                         5d9s22
      dimension i4o(*),ionex(*),noc(*),nvirt(*),isblk(4,*),isblkx(4,*)  5d9s22
      include "common.store"
      logical dob                                                       7d18s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
c
c     derivative integrals were computed full, replicated across all
c     procs. Now
c     a) only store this procs columns, and
c     b) only store row triangle
c
      do is=1,nsblk
       call ilimts(noc(isblk(3,is)),noc(isblk(4,is)),mynprocg,          5d31s22
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          5d9s22
       nhere=ih+1-il                                                    5d9s22
       if(nhere.gt.0)then                                               5d9s22
c
c     a)
c
        nrowo=noc(isblk(1,is))*noc(isblk(2,is))                         5d31s22
        ncol=noc(isblk(3,is))*noc(isblk(4,is))
        if(idwsdeb.gt.100)then                                          5d31s22
         write(6,*)('for '),(isblk(j,is),j=1,4)
         write(6,*)('go from ')
         call prntm2(bc(i4o(is)),nrowo,ncol,nrowo)
        end if
        ngot=nrowo*nhere                                                5d31s22
        istrt=i4o(is)+nrowo*(il-1)                                      5d31s22
        ngotm=ngot-1                                                    5d31s22
        do i=0,ngotm                                                    5d31s22
         bc(i4o(is)+i)=bc(istrt+i)                                       5d31s22
        end do                                                          5d31s22
        if(idwsdeb.gt.100)then                                          5d31s22
         write(6,*)('to ')
         call prntm2(bc(i4o(is)),nrowo,nhere,nrowo)
        end if                                                          5d31s22
c
c     b)
c
        if(isblk(1,is).eq.isblk(2,is).and.dob)then                      7d18s22
         if(min(noc(isblk(1,is)),noc(isblk(3,is))).gt.0)then             5d9s22
          nrown=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                5d9s22
          if(idwsdeb.gt.10)then                                         5d19s22
           write(6,*)('triangularizing 4o type '),(isblk(j,is),j=1,4)
           write(6,*)('go from ')
           call prntm2(bc(i4o(is)),nrowo,nhere,nrowo)
          end if                                                        5d19s22
          itmp=ibcoff                                                   5d9s22
          ibcoff=itmp+nhere*nrown                                       5d9s22
          call enough('trianglez.  1',bc,ibc)
          jtmp=itmp                                                     5d9s22
          do icol=0,nhere-1
           do i2=0,noc(isblk(1,is))-1                                   5d9s22
            i12=i4o(is)+noc(isblk(1,is))*(i2+noc(isblk(1,is))*icol)     5d9s22
            do i1=0,i2                                                  5d9s22
             bc(jtmp+i1)=bc(i12+i1)                                     5d9s22
            end do                                                      5d9s22
            jtmp=jtmp+i2+1                                              5d9s22
           end do                                                       5d9s22
          end do
          do i=0,nrown*nhere-1                                          5d9s22
           bc(i4o(is)+i)=bc(itmp+i)                                     5d9s22
          end do                                                        5d9s22
          if(idwsdeb.gt.10)then                                         5d19s22
           write(6,*)('to ')
           call prntm2(bc(i4o(is)),nrown,nhere,nrown)                    5d9s22
          end if                                                        5d19s22
          ibcoff=itmp                                                   5d9s22
         end if                                                         5d9s22
        end if                                                          5d9s22
       end if                                                           5d9s22
      end do
      do is=1,nsblkx
       nrowo=noc(isblkx(1,is))*noc(isblkx(2,is))                        5d31s22
       call ilimts(noc(isblkx(3,is)),nvirt(isblkx(4,is)),mynprocg,      5d9s22
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          5d9s22
       nhere=ih+1-il                                                    5d9s22
       if(nhere.gt.0)then                                               5d31s22
c
c     a)                                                                5d31s22
c                                                                       5d31s22
        istrt=ionex(is)+nrowo*(il-1)                                    5d31s22
        ngot=nrowo*nhere                                                5d31s22
        ngotm=ngot-1                                                    5d31s22
        do i=0,ngotm                                                    5d31s22
         bc(ionex(is)+i)=bc(istrt+i)                                    5d31s22
        end do                                                          5d31s22
        if(isblkx(1,is).eq.isblkx(2,is).and.dob)then                    7d18s22
c
c     b)
c
         if(min(noc(isblkx(1,is)),noc(isblkx(3,is))).gt.0)then             5d9s22
          if(nvirt(isblkx(4,is)).gt.0)then                                5d9s22
           nrown=(noc(isblkx(1,is))*(noc(isblkx(1,is))+1))/2                5d9s22
           if(idwsdeb.gt.10)then                                        5d19s22
            write(6,*)('triangularizing onex type '),
     $          (isblkx(j,is),j=1,4)
            write(6,*)('go from ')
            call prntm2(bc(ionex(is)),nrowo,nhere,nrowo)
           end if                                                       5d19s22
           itmp=ibcoff                                                   5d9s22
           ibcoff=itmp+nhere*nrown                                       5d9s22
           call enough('trianglez.  2',bc,ibc)
           jtmp=itmp                                                     5d9s22
           do icol=0,nhere-1
            do i2=0,noc(isblkx(1,is))-1                                   5d9s22
             i12=ionex(is)+noc(isblkx(1,is))*(i2+noc(isblkx(1,is))*icol)  5d9s22
             do i1=0,i2                                                  5d9s22
              bc(jtmp+i1)=bc(i12+i1)                                     5d9s22
             end do                                                      5d9s22
             jtmp=jtmp+i2+1                                              5d9s22
            end do                                                       5d9s22
           end do
           do i=0,nrown*nhere-1                                          5d9s22
            bc(ionex(is)+i)=bc(itmp+i)                                     5d9s22
           end do                                                        5d9s22
           if(idwsdeb.gt.10)then                                        5d19s22
            write(6,*)('to ')
            call prntm2(bc(ionex(is)),nrown,nhere,nrown)                    5d9s22
           end if                                                       5d19s22
           ibcoff=itmp                                                   5d9s22
          end if                                                         5d9s22
         end if                                                         5d9s22
        end if                                                          5d9s22
       end if                                                           5d9s22
      end do
      return
      end
