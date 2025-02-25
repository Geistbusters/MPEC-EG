c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine taylor(amat,umat,no,ni,iret,bc,ibc)                    11d9s22
      implicit real*8 (a-h,o-z)
      dimension amat(no,ni),umat(*)                                     1d3s18
      include "common.store"                                            1d3s18
      iret=0                                                            1d3s18
      nb=no+ni                                                          1d3s18
      thres=1d-10
      it1=ibcoff                                                        1d3s18
      it2=it1+nb*nb                                                     1d3s18
      it3=it2+nb*nb                                                     1d3s18
      it4=it3+nb*nb                                                     1d3s18
      nbm=nb-1                                                          1d3s18
      ibcoff=it4+nb*nb                                                  1d3s18
      call enough('taylor.  1',bc,ibc)
      jt3=it3                                                           1d3s18
      jt1=it1                                                           1d3s18
      do i=0,nbm                                                        1d3s18
       do j=0,nbm                                                       1d3s18
        bc(jt3+j)=0d0                                                   1d3s18
        bc(jt1+j)=0d0                                                   1d3s18
       end do                                                           1d3s18
       bc(jt3+i)=1d0                                                    1d3s18
       jt3=jt3+nb                                                       1d3s18
       jt1=jt1+nb                                                       1d3s18
      end do                                                            1d3s18
      do i=1,ni                                                         1d3s18
       im=i-1                                                           1d3s18
       do j=1,no                                                        1d3s18
        jp=j+ni-1                                                       1d3s18
        ji=it3+jp+nb*im                                                 1d3s18
        ij=it3+im+nb*jp                                                 1d3s18
        bc(ji)=amat(j,i)                                                1d3s18
        bc(ij)=-amat(j,i)                                               1d3s18
        ji=it1+jp+nb*im                                                 1d3s18
        ij=it1+im+nb*jp                                                 1d3s18
        bc(ji)=amat(j,i)                                                1d3s18
        bc(ij)=-amat(j,i)                                               1d3s18
       end do                                                           1d3s18
      end do                                                            1d3s18
      nbbm=nb*nb-1                                                      1d3s18
      do i=0,nbbm                                                       1d3s18
       bc(it2+i)=bc(it1+i)                                              1d3s18
      end do                                                            1d3s18
      fact=0.5d0                                                        1d3s18
      fact2=2d0                                                         1d3s18
      iter=0                                                            1d3s18
    1 continue                                                          1d3s18
       iter=iter+1                                                      1d3s18
       call dgemm('n','n',nb,nb,nb,1d0,bc(it1),nb,bc(it2),nb,0d0,       1d3s18
     $      bc(it4),nb,                                                 1d3s18
     d' taylor.  1')
       do i=0,nbbm                                                      1d3s18
        bc(it2+i)=bc(it4+i)                                             1d3s18
        bc(it3+i)=bc(it3+i)+bc(it2+i)*fact                              1d3s18
       end do                                                           1d3s18
       do i=0,nbm                                                       1d3s18
        bc(it4+i)=0d0                                                   1d3s18
       end do                                                           1d3s18
       jt3=it3                                                          1d3s18
       do i=0,nbm                                                       1d3s18
        do j=0,nbm                                                      1d3s18
         bc(it4+j)=bc(it4+j)+bc(jt3+j)**2                               1d3s18
        end do                                                          1d3s18
        jt3=jt3+nb                                                      1d3s18
       end do                                                           1d3s18
       rms=0d0                                                          1d3s18
       do i=0,nbm                                                       1d3s18
        rms=rms+(bc(it4+i)-1d0)**2                                      1d3s18
       end do                                                           1d3s18
       rms=sqrt(rms/dfloat(nb))                                         1d3s18
       if(rms.gt.thres.and.iter.lt.15)then                               5d22s23
        fact2=fact2+1d0                                                  1d3s18
        fact=fact/fact2                                                  1d3s18
        go to 1                                                          1d3s18
       end if                                                           1d3s18
      call dgemm('n','n',nb,nb,nb,1d0,umat,nb,bc(it3),nb,0d0,bc(it1),nb,1d3s18
     d' taylor.  2')
      do i=0,nbm                                                        1d3s18
       iiv=it1+nb*i                                                     1d3s18
       do j=0,i-1                                                       1d3s18
        jjv=it1+nb*j                                                    1d3s18
        dot=0d0                                                         1d3s18
        do k=0,nbm                                                      1d3s18
         dot=dot+bc(iiv+k)*bc(jjv+k)                                    1d3s18
        end do                                                          1d3s18
        do k=0,nbm                                                      1d3s18
         bc(iiv+k)=bc(iiv+k)-dot*bc(jjv+k)                              1d3s18
        end do                                                          1d3s18
       end do                                                           1d3s18
       sum=0d0                                                          1d3s18
       do k=0,nbm                                                       1d3s18
        sum=sum+bc(iiv+k)**2                                            1d3s18
       end do                                                           1d3s18
       sum=1d0/sqrt(sum)                                                1d3s18
       do k=0,nbm                                                       1d3s18
        bc(iiv+k)=bc(iiv+k)*sum                                         1d3s18
       end do                                                           1d3s18
      end do                                                            1d3s18
      do i=0,nbbm                                                       1d3s18
       ip=i+1                                                           1d3s18
       umat(ip)=bc(it1+i)                                               1d3s18
      end do                                                            1d3s18
      ibcoff=it1                                                        1d3s18
      return                                                            1d3s18
      end                                                               1d3s18
