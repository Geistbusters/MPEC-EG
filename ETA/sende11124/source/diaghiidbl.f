c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diaghiidbl(iaorb,nalpha,numa,iborb,nbeta,numb,         7d21s22
     $     jdenpt,nsymb,nsbeta,ivec,ilc,ihc,bc,ibc,nroot,mdenoff,       3d15s23
     $     iden1,igoal)                                                       3d15s23
      implicit real*8 (a-h,o-z)
c
c     compute bl density from diagonal matrix elements of hamiltonian
c
      include "common.store"
      include "common.cas"
      integer*1 iaorb(nalpha,numa),iborb(nbeta,numb)
      dimension ivec(*),jdenpt(*),ilc(*),ihc(*),nsbeta(*),iden1(*)      3d15s23
      norb=iacto(1)                                                     8d8s06
      do i=2,nsymb                                                      8d8s06
       norb=norb+iacto(i)                                               8d8s06
      end do                                                            8d8s06
      isma=ibcoff                                                        8d8s06
      ireloa=isma+norb                                                  8d8s06
      ivv=ireloa+norb                                                   3d14s23
      ntri=(nroot*(nroot+1))/2                                          3d14s23
      ibcoff=ivv+ntri                                                   3d14s23
      call enough('diaghiidbl.  1',bc,ibc)
      jsm=isma-1                                                         8d8s06
      jrelo=ireloa-1                                                     8d8s06
      do i=1,nsymb                                                      8d8s06
       do j=1,iacto(i)                                                  8d8s06
        ibc(jsm+j)=i                                                    8d8s06
        ibc(jrelo+j)=j                                                  8d8s06
  153   format(3i5)
       end do
       jsm=jsm+iacto(i)                                                 8d8s06
       jrelo=jrelo+iacto(i)                                             8d8s06
      end do                                                            8d8s06
      ismb=ibcoff                                                        8d8s06
      irelob=ismb+norb                                                  8d8s06
      ibcoff=irelob+norb                                                8d8s06
      call enough('diaghiidbl.  2',bc,ibc)
      jsm=ismb-1                                                         8d8s06
      jrelo=irelob-1                                                     8d8s06
      do i=1,nsymb                                                      8d8s06
       do j=1,iacto(i)                                                  8d8s06
        ibc(jsm+j)=i                                                    8d8s06
        ibc(jrelo+j)=j                                                  8d8s06
       end do
       jsm=jsm+iacto(i)                                                 8d8s06
       jrelo=jrelo+iacto(i)                                             8d8s06
      end do                                                            8d8s06
c
c     alpha-beta contributions
c
      ii=0
      jreloa=ireloa-1                                                   8d8s06
      jrelob=irelob-1                                                   8d8s06
      jsma=isma-1                                                       8d8s06
      jsmb=ismb-1                                                       8d8s06
      ioffa=0                                                           8d8s06
      do isb=1,nsymb                                                    8d8s06
       ioffb=0                                                          8d8s06
       nha=ihc(isb)+1-ilc(isb)                                          7d25s22
       do is2=1,nsbeta(isb)-1                                           8d8s06
        ioffb=ioffb+nbdet(is2)                                          8d8s06
       end do                                                           8d8s06
       do ia=ilc(isb),ihc(isb)                                          8d8s06
        iap=ia+ioffa                                                    8d8s06
        do ib=1,nbdet(nsbeta(isb))                                      8d8s06
         ibp=ib+ioffb                                                   8d8s06
         ii=ii+1
         iad=ivec(isb)+nroot*(ia-ilc(isb)+nha*(ib-1))                   3d14s23
         jvv=ivv                                                        3d14s23
         do k1=0,nroot-1                                                3d14s23
          do k2=0,k1                                                    3d14s23
           bc(jvv)=-2d0*bc(iad+k1)*bc(iad+k2)                           3d14s23
           jvv=jvv+1                                                    3d14s23
          end do                                                        3d14s23
         end do                                                         3d14s23
c
c     for diagonals, we have sum a<b -2*abba
c
         do j=1,nalpha
          is14=ibc(jsma+iaorb(j,iap))                                              8d8s06
          if14=ibc(jreloa+iaorb(j,iap))-1                               7d21s22
          do k=0,ntri-1                                                 3d15s23
           kp=k+mdenoff                                                 3d15s23
           iad=iden1(is14)+if14+iacto(is14)*(if14+iacto(is14)*kp)       3d15s23
           bc(iad)=bc(iad)-0.5d0*bc(ivv+k)                              3d15s23
          end do                                                        3d15s23
          do jp=j+1,nalpha                                              7d21s22
           is23=ibc(jsma+iaorb(jp,iap))                                           8d8s06
           if23=ibc(jreloa+iaorb(jp,iap))-1                             7d21s22
           i2eu=ifind2(is14,is23,is23,is14,icase)
           if(icase.eq.1)then                                             7d7s22
            do k=0,ntri-1                                               3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+if14+iacto(is14)*(if23+iacto(is23)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                       3d14s23
             bc(iad)=bc(iad)+bc(ivv+k)                                  3d14s23
            end do                                                      3d14s23
           else if(icase.eq.4)then                                        7d7s22
            do k=0,ntri-1                                               3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+if23+iacto(is23)*(if14+iacto(is14)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                       3d14s23
             bc(iad)=bc(iad)-bc(ivv+k)                                  3d14s23
            end do                                                      3d14s23
           else if(icase.eq.3)then                                        7d7s22
            do k=0,ntri-1                                               3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+if14+iacto(is14)*(if23+iacto(is23)*(if23  3d14s23
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                       3d14s23
             bc(iad)=bc(iad)-bc(ivv+k)                                  3d14s23
            end do                                                      3d14s23
           else if(icase.eq.2)then                                        7d7s22
            do k=0,ntri-1                                               3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+if23+iacto(is23)*(if14+iacto(is14)*(if14     7d7s22
     $        +iacto(is14)*(if23+iacto(is23)*kp)))                       3d14s23
             bc(iad)=bc(iad)+bc(ivv+k)                                  3d14s23
            end do                                                      3d14s23
           end if                                                         7d7s22
          end do
         end do
         do j=1,nbeta
          is14=ibc(jsmb+iborb(j,ibp))                                              8d8s06
          if14=ibc(jrelob+iborb(j,ibp))-1                               7d21s22
          do k=0,ntri-1                                                 3d15s23
           kp=k+mdenoff                                                 3d15s23
           iad=iden1(is14)+if14+iacto(is14)*(if14+iacto(is14)*kp)       3d15s23
           bc(iad)=bc(iad)-0.5d0*bc(ivv+k)                              3d15s23
          end do                                                        3d15s23
          do jp=j+1,nbeta                                               7d21s22
           is23=ibc(jsmb+iborb(jp,ibp))                                           8d8s06
           if23=ibc(jrelob+iborb(jp,ibp))-1                             7d21s22
           i2eu=ifind2(is14,is23,is23,is14,icase)
           if(icase.eq.1)then                                             7d7s22
            do k=0,ntri-1                                               3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+if14+iacto(is14)*(if23+iacto(is23)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                       3d14s23
             bc(iad)=bc(iad)+bc(ivv+k)                                  3d14s23
            end do                                                      3d14s23
           else if(icase.eq.4)then                                        7d7s22
            do k=0,ntri-1                                               3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+if23+iacto(is23)*(if14+iacto(is14)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                       3d14s23
             bc(iad)=bc(iad)-bc(ivv+k)                                  3d14s23
            end do                                                      3d14s23
           else if(icase.eq.3)then                                        7d7s22
            do k=0,ntri-1                                               3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+if14+iacto(is14)*(if23+iacto(is23)*(if23     7d7s22
     $        +iacto(is23)*(if14+iacto(is14)*kp)))                       3d14s23
             bc(iad)=bc(iad)-bc(ivv+k)                                  3d14s23
            end do                                                      3d14s23
           else if(icase.eq.2)then                                        7d7s22
            do k=0,ntri-1                                               3d14s23
             kp=k+mdenoff                                                 3d14s23
             iad=jdenpt(i2eu)+if23+iacto(is23)*(if14+iacto(is14)*(if14     7d7s22
     $        +iacto(is14)*(if23+iacto(is23)*kp)))                       3d14s23
             bc(iad)=bc(iad)+bc(ivv+k)                                  3d14s23
            end do                                                      3d14s23
           end if                                                         7d7s22
          end do
         end do
        end do
       end do
       ioffa=ioffa+nadet(isb)                                           8d8s06
      end do                                                            8d8s06
      ibcoff=isma                                                       8d8s06
      return
      end
