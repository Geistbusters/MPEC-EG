c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hessycmo(nv4,idav,ibig,iyc,noca,nocb,nvirt,            9d16s04
     $     nvirtb,ibmat1,ibmat2,igmat,i2s,i2e,i1s,i1e,xlmi2,iflg,bc,ibc)11d10s22
      implicit real*8 (a-h,o-z)
      integer iarg1,iarg2,iarg3,nproc,idum,mtmp,iarg4                   3d12s12
      include "common.store"
      jbig=ibig+nv4*(idav-1)
      bc(jbig)=0d0
      iaohv=ibcoff
      ibcoff=iaohv+nocb*nvirtb                                          9d16s04
      iaohvt=ibcoff
      ibcoff=iaohvt+nocb*nvirtb                                         9d16s04
      call enough('hessycmo.  1',bc,ibc)
      jbmat2=ibmat2-1
      fact=2d0*bc(iyc)
      jaohv=iaohv-1
      jaohvt=iaohvt-1
       if(nv4.ne.0)then                                                 12d6s05
        do i=1,nvirtb*nocb                                              12d6s05
         bc(jbig)=bc(jbig)+2d0*bc(iyc+i)*bc(jbmat2+i)
         bc(jaohv+i)=fact*bc(jbmat2+i)
         bc(jaohvt+i)=0d0
        end do                                                          12d6s05
       else                                                             12d6s05
        do i=1,nvirtb*nocb                                              12d6s05
         bc(jaohvt+i)=0d0                                               12d6s05
         bc(jaohv+i)=0d0                                                12d6s05
        end do                                                          12d6s05
       end if                                                           12d6s05
       if(iflg.gt.0)then                                                3d12s12
        ii=0
        jgmat=igmat
        do i2=i2s,i2e
         if(i2.eq.i2s)then
          i10=i1s
         else
          i10=1
         end if
         if(i2.eq.i2e)then
          i1n=i1e
         else
          i1n=nvirtb
         end if
         do i1=i10,i1n                                                   3d16s07
          jaov=iyc+i2-nvirt                                              3d16s07
          do i3=1,nocb                                                   7d26s07
           jaohvt=iaohvt+i3-1+(i1-1)*nocb                                7d26s07
           do i4=1,noca                                                  7d26s07
            fact=2d0*xlmi2*bc(jaov+i4*nvirt)                             7d26s07
            bc(jaohvt)=bc(jaohvt)+fact*bc(jgmat)                         3d16s07
            jgmat=jgmat+1
           end do                                                        3d16s07
          end do                                                         3d16s07
          if(i1.eq.i2.and.nv4.ne.0)then
           do i3=1,nocb
            jaohvt=iaohvt+i3-1+(i1-1)*nocb                                3d16s07
            j1=ibmat1-1+(i3-1)*nocb
            j2=ibmat1+i3-1-nocb
            do i4=1,nocb
             fact=bc(jaov+i4*nvirtb)*xlmi2
             bc(jaohvt)=bc(jaohvt)
     $            -(bc(j1+i4)+bc(j2+i4*nocb))*fact
            end do                                                       3d16s07
           end do                                                        3d16s07
          end if
         end do                                                          3d16s07
        end do
       else                                                             3d12s12
        ii=0
        jgmat=igmat
        do i2=i2s,i2e
         if(i2.eq.i2s)then
          i10=i1s
         else
          i10=1
         end if
         if(i2.eq.i2e)then
          i1n=i1e
         else
          i1n=nvirt                                                     3d12s12
         end if
         do i1=i10,i1n                                                   3d16s07
          jaov=iyc+i1-nvirt                                             3d12s12
          do i4=1,noca                                                  3d12s12
           fact=2d0*xlmi2*bc(jaov+i4*nvirt)                              3d13s12
           do i3=1,nocb                                                 3d12s12
            jaohvt=iaohvt+i3-1+(i2-1)*nocb                              3d13s12
            bc(jaohvt)=bc(jaohvt)+fact*bc(jgmat)                         3d16s07
            jgmat=jgmat+1
           end do                                                        3d16s07
          end do                                                         3d16s07
          if(i1.eq.i2.and.nv4.ne.0)then
           do i3=1,nocb
            jaohvt=iaohvt+i3-1+(i1-1)*nocb                                3d16s07
            j1=ibmat1-1+(i3-1)*nocb
            j2=ibmat1+i3-1-nocb
            do i4=1,nocb
             fact=bc(jaov+i4*nvirt)*xlmi2                               3d13s12
             bc(jaohvt)=bc(jaohvt)
     $            -(bc(j1+i4)+bc(j2+i4*nocb))*fact
            end do                                                       3d16s07
           end do                                                        3d16s07
          end if
         end do                                                          3d16s07
        end do
       end if                                                           3d12s12
       iarg1=nocb*nvirtb                                                3d16s07
       call dws_gsumf(bc(iaohvt),iarg1)                                 3d12s12
      do 10 i=1,nocb
       jaohv=iaohv-1+(i-1)*nvirtb
       jaohvt=iaohvt+i-1-nocb
       kbig=jbig+(i-1)*nvirtb
       do 11 j=1,nvirtb
        bc(kbig+j)=bc(jaohv+j)+bc(jaohvt+j*nocb)
   11  continue
   10 continue
      ibcoff=iaohv
      return
      end
