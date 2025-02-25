c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine diaghiid(vec1,vec2,nconf,iaorb,nalpha,numa,            5d16s22
     $                   iborb,nbeta,numb,iden1,jdenpt,shift,ilc,ihc,       8d8s06
     $                   debug,nsymb,nsbeta,nroot,wgt,l2e,bc,ibc,       3d15s23
     $     mdenoff,igoal)                                                     3d15s23
      implicit real*8 (a-h,o-z)
c
c     compute density from diagonal matrix elements of hamiltonian
c
      include "common.store"
      include "common.cas"
      integer*1 iaorb(nalpha,numa),iborb(nbeta,numb)
      integer*8 iarg1,iarg2                                             5d7s18
      integer*2 ihdig(4,nconf)                                          8d22s06
      logical l2e                                                       6d2s22
      dimension vec1(nconf,nroot),jdenpt(1),ilc(8),ihc(8),nsbeta(8),    5d16s22
     $     iden1(8),wgt(nroot),vec2(nconf,nroot)                        5d16s22
      norb=iacto(1)                                                     8d8s06
      do i=2,nsymb                                                      8d8s06
       norb=norb+iacto(i)                                               8d8s06
      end do                                                            8d8s06
      isma=ibcoff                                                        8d8s06
      ireloa=isma+norb                                                  8d8s06
      ibcoff=ireloa+norb                                                8d8s06
      call enough('diaghiid.  1',bc,ibc)
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
      call enough('diaghiid.  2',bc,ibc)
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
      ivv=ibcoff                                                        3d15s23
      if(.not.l2e)then                                                  3d15s23
       ntri=(nroot*(nroot+1))/2                                         3d15s23
       ibcoff=ivv+ntri                                                  3d15s23
       call enough('diaghiid.ivv',bc,ibc)                               3d15s23
      end if                                                            3d15s23
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
       do is2=1,nsbeta(isb)-1                                           8d8s06
        ioffb=ioffb+nbdet(is2)                                          8d8s06
       end do                                                           8d8s06
       do ia=ilc(isb),ihc(isb)                                          8d8s06
        iap=ia+ioffa                                                    8d8s06
        do ib=1,nbdet(nsbeta(isb))                                      8d8s06
         ibp=ib+ioffb                                                   8d8s06
         ii=ii+1
         if(l2e)then                                                    3d15s23
          vv=0d0                                                         8d19s14
          do ir=1,nroot                                                  8d19s14
           vv=vv+wgt(ir)*vec1(ii,ir)*vec2(ii,ir)                         5d16s22
          end do                                                         8d19s14
         else                                                           3d15s23
          jvv=ivv                                                       3d15s23
          do k1=1,nroot                                                 3d15s23
           do k2=1,k1                                                   3d15s23
            bc(jvv)=vec1(ii,k1)*vec2(ii,k2)                             3d15s23
            jvv=jvv+1                                                   3d15s23
           end do                                                       3d15s23
          end do                                                        3d15s23
         end if                                                         3d15s23
         do j=1,nalpha
          jo=ibc(jreloa+iaorb(j,iap))                                            8d8s06
          jj=((jo*(jo+1))/2)-1
          nj=ibc(jsma+iaorb(j,iap))                                              8d8s06
          na=(iacto(nj)*(iacto(nj)+1))/2                                  8d8s06
          if(l2e)then                                                   3d15s23
           iadd=iden1(nj)+(jo-1)*(iacto(nj)+1)                              8d8s06
           bc(iadd)=bc(iadd)+vv
          else                                                          3d15s23
           iadd=iden1(nj)+jo-1+iacto(nj)*(jo-1+iacto(nj)*mdenoff)       3d15s23
           nadd=iacto(nj)*iacto(nj)                                     3d15s23
           do k1=1,nroot                                                3d27s23
            do k2=1,k1                                                  3d27s23
             term=vec1(ii,k1)*vec2(ii,k2)                               3d27s23
             bc(iadd)=bc(iadd)+term
             iadd=iadd+nadd                                              3d15s23
            end do                                                      3d27s23
           end do                                                       3d27s23
          end if                                                        3d15s23
          if(l2e)then                                                   6d2s22
           do jp=1,nalpha
            jpo=ibc(jreloa+iaorb(jp,iap))                                         8d8s06
            jjp=((jpo*(jpo+1))/2)-1
            njp=ibc(jsma+iaorb(jp,iap))                                           8d8s06
            i2eu=ifind2(nj,nj,njp,njp,icase)                               8d15s06
            iadd1=jdenpt(i2eu)+jj+jjp*na                                      8d8s06
            bc(iadd1)=bc(iadd1)+0.5d0*vv                                 8d19s14
            i2eu=ifind2(nj,njp,nj,njp,icase)                               8d15s06
            if(icase.eq.2)then                                             8d15s06
             iadd2=jdenpt(i2eu)+jpo-1+iacto(njp)*(jo-1+iacto(nj)              8d15s06
     $         *(jpo-1+iacto(njp)*(jo-1)))
            else if(nj.eq.njp)then
             ix=max(jpo,jo)
             in=min(jpo,jo)
             kk=((ix*(ix-1))/2)+in-1
             iadd2=jdenpt(i2eu)+kk*(na+1)                                     8d14s06
            else                                                           8d8s06
             iadd2=jdenpt(i2eu)+jo-1+iacto(nj)*(jpo-1+iacto(njp)              8d14s06
     $           *(jo-1+iacto(nj)*(jpo-1)))                               8d8s06
            end if                                                         8d8s06
            bc(iadd2)=bc(iadd2)-0.5d0*vv
           end do
          end if                                                        6d2s22
         end do
         do j=1,nbeta
          jo=ibc(jrelob+iborb(j,ibp))                                            8d8s06
          jj=((jo*(jo+1))/2)-1
          nj=ibc(jsmb+iborb(j,ibp))                                              8d8s06
          na=(iacto(nj)*(iacto(nj)+1))/2                                  8d8s06
          if(l2e)then                                                   3d15s23
           iadd=iden1(nj)+(jo-1)*(iacto(nj)+1)                              8d8s06
           bc(iadd)=bc(iadd)+vv
          else                                                          3d15s23
           iadd=iden1(nj)+jo-1+iacto(nj)*(jo-1+iacto(nj)*mdenoff)       3d15s23
           nadd=iacto(nj)*iacto(nj)                                     3d15s23
           do k1=1,nroot                                                3d27s23
            do k2=1,k1                                                  3d27s23
             term=vec1(ii,k1)*vec2(ii,k2)                               3d27s23
             bc(iadd)=bc(iadd)+term
             iadd=iadd+nadd                                              3d15s23
            end do                                                      3d27s23
           end do                                                       3d27s23
          end if                                                        3d15s23
          if(l2e)then                                                   6d2s22
           do jp=1,nbeta
            jpo=ibc(jrelob+iborb(jp,ibp))                                         8d8s06
            jjp=((jpo*(jpo+1))/2)-1
            njp=ibc(jsmb+iborb(jp,ibp))                                           8d8s06
            i2eu=ifind2(nj,nj,njp,njp,icase)                               8d15s06
            iadd1=jdenpt(i2eu)+jj+jjp*na                                      8d8s06
            bc(iadd1)=bc(iadd1)+0.5d0*vv                                 8d19s14
            i2eu=ifind2(nj,njp,nj,njp,icase)                               8d15s06
            if(icase.eq.2)then                                             8d15s06
             iadd2=jdenpt(i2eu)+jpo-1+iacto(njp)*(jo-1+iacto(nj)              8d15s06
     $         *(jpo-1+iacto(njp)*(jo-1)))
            else if(nj.eq.njp)then
             ix=max(jpo,jo)
             in=min(jpo,jo)
             kk=((ix*(ix-1))/2)+in-1
             iadd2=jdenpt(i2eu)+kk*(na+1)                                     8d14s06
            else                                                           8d8s06
             iadd2=jdenpt(i2eu)+jo-1+iacto(nj)*(jpo-1+iacto(njp)              8d14s06
     $           *(jo-1+iacto(nj)*(jpo-1)))                               8d8s06
            end if                                                         8d8s06
            bc(iadd2)=bc(iadd2)-0.5d0*vv
           end do
          end if                                                        6d2s22
         end do
         if(l2e)then                                                    6d2s22
          do ja=1,nalpha
           jao=ibc(jreloa+iaorb(ja,iap))                                 8d8s06
           jja=((jao*(jao+1))/2)-1
           na=ibc(jsma+iaorb(ja,iap))                                    8d8s06
           nna=(iacto(na)*(iacto(na)+1))/2                               8d8s06
           do jb=1,nbeta
            jbo=ibc(jrelob+iborb(jb,ibp))                                8d8s06
            jjb=((jbo*(jbo+1))/2)-1
            nb=ibc(jsmb+iborb(jb,ibp))                                   8d8s06
            i2eu=ifind2(nb,nb,na,na,icase)                               8d15s06
            if(icase.eq.1)then                                           4d24s07
             nna=(iacto(nb)*(iacto(nb)+1))/2                              8d14s06
             iadd=jdenpt(i2eu)+jjb+jja*nna                                  4d24s07
            else if(icase.eq.5)then                                      4d24s07
             nna=(iacto(na)*(iacto(na)+1))/2                             4d24s07
             iadd=jdenpt(i2eu)+jja+jjb*nna                                  4d24s07
            else                                                         4d24s07
             write(6,*)('icase is not 1 or 5 '),icase,nb,na
             stop
            end if
            iadd=jdenpt(i2eu)+jjb+jja*nna
            bc(iadd)=bc(iadd)+vv                                         8d19s14
           end do
          end do
         end if                                                         6d2s22
 9332    format(i4,2i5,i2,3f18.8)
        end do
       end do
       ioffa=ioffa+nadet(isb)                                           8d8s06
      end do                                                            8d8s06
      ibcoff=isma                                                       8d8s06
      return
      end
