c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine gotoqs(va,vq,ndima,ndimq,nroot,ipoint,myfcn,ibas,      3d25s21
     $     ncsf,mdon)                                                   3d25s21
      implicit real*8 (a-h,o-z)                                         3d25s21
      dimension va(ndima,*),vq(ndimq,*),ipoint(*),ibas(3,*),ncsf(*)     3d25s21
      ioffa=0                                                           3d25s21
      ioffq=0                                                           3d25s21
      do if=1,myfcn                                                     3d25s21
       nclo=ibas(1,if)                                                  3d25s21
       nclop=nclo+1                                                     3d25s21
       iarg=nclop-mdon                                                  3d25s21
       if(ipoint(if).eq.1)then                                          3d25s21
        do ir=1,nroot                                                   3d25s21
         do i=1,ncsf(iarg)                                              3d25s21
          vq(ioffq+i,ir)=va(ioffa+i,ir)                                 3d25s21
         end do                                                         3d25s21
        end do                                                          3d25s21
        ioffq=ioffq+ncsf(iarg)                                          3d25s21
       end if                                                           3d25s21
       ioffa=ioffa+ncsf(iarg)                                           3d25s21
      end do                                                            3d25s21
      return                                                            3d25s21
      end                                                               3d25s21
