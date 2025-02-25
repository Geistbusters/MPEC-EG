c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine ssdvrc(lb,b1,b2,b3,b4,b5,lk,k1,k2,k3,k4,k5,iout,data,   3d4s20
     $     iflag,clight,bc,ibc)                                         11d9s22
      implicit real*8 (a-h,o-z)
      real*8 k1,k2,k3,k4,k5
c
c     call either onep
c     or derid to generate 2 component matrix element for               3d4s20
c     spin-free basis in full 4 component space.
c     if iflag=0, we are doing small-small overlap, data not used.      3d4s20
c     if iflag=2, we are doing small-small nuclear attraction,          3d4s20
c     data contains nuclear position and size.                          3d4s20
c
      dimension data(*),ider(3,3),nspin(2,2),ispin(3,2,2,2),phs(2),     3d4s20
     $     ioffb(2),ioffk(2),cvec(6,2,2),ip2eri(3,3)                    3d25s20
      data ider/1,0,0, 0,1,0, 0,0,1/                                    3d4s20
      data phs/1d0,-1d0/                                                3d4s20
      data nspin/1,2,2,1/
      include "common.store"                                            3d4s20
      ibcoffo=ibcoff                                                    3d4s20
      nlb=2*lb+1                                                        3d4s20
      nlb2=nlb*2                                                        3d4s20
      nlk=2*lk+1                                                        3d4s20
      nlk2=nlk*2                                                        3d4s20
      iout=ibcoff                                                       3d4s20
      iouti=iout+nlb2*nlk2                                              3d4s20
      ibcoff=iouti+nlb2*nlk2                                            3d4s20
      iprim=ibcoff                                                      3d25s20
      ibcoff=iprim+nlb2*nlk2*9                                          3d25s20
      call enough('ssdvrc.  1',bc,ibc)
      do i=iout,ibcoff-1                                                3d4s20
       bc(i)=0d0                                                        3d4s20
      end do                                                            3d4s20
      fact=0.5d0/clight                                                 3d25s20
c
c     contraction from real 4 component basis functions to complex
c     2 component basis functions
c
      do i3=1,2                                                         1d13s23
       do i2=1,2                                                        1d13s23
        do i1=1,2                                                       1d13s23
         cvec(i1,i2,i3)=0d0                                             1d13s23
        end do                                                          1d13s23
       end do                                                           1d13s23
      end do                                                            1d13s23
c     what dyall used
      cvec(3,1,2)=-fact
      cvec(4,1,2)=-fact
      cvec(5,1,1)=fact
      cvec(1,2,2)=-fact
      cvec(2,2,1)=-fact
      cvec(6,2,2)=fact
      do ik=1,3                                                         3d25s20
       do ib=1,3                                                        3d25s20
        if(iflag.eq.0)then                                              3d25s20
         call onep(lb,b1,b2,b3,b4,b5,lk,k1,k2,k3,k4,k5,idum,ibo,        3d25s20
     $           ider(1,ib),ider(2,ib),ider(3,ib),0,0,0,                3d25s20
     $           ider(1,ik),ider(2,ik),ider(3,ik),bc,ibc)               11d9s22
        else                                                            3d25s20
         call derid(lb,b1,b2,b3,b4,b5,idum,lk,k1,k2,k3,k4,k5,           3d25s20
     $           idum,0,data(1),data(2),data(3),data(4),data(5),idum,   3d25s20
     $           0,data(1),data(2),data(3),data(4),data(5),idum,dum,1,  3d25s20
     $           idum,ibo,ider(1,ib),ider(2,ib),ider(3,ib),             3d25s20
     $           ider(1,ik),ider(2,ik),ider(3,ik),0,0,0, 0,0,0, 1,      3d25s20
     $           .false.,data(6),bc,ibc)                                11d14s22
        end if                                                          3d25s20
        ip2eri(ib,ik)=ibo                                               3d25s20
        ibcoff=ibo+nlb*nlk                                              3d25s20
       end do                                                           3d25s20
      end do                                                            3d25s20
      nlk3=nlk*3                                                        3d25s20
      nlb3=nlb*3                                                        3d25s20
      nlb6=nlb3*2                                                       3d25s20
      do ispinks=0,1                                                    3d25s20
       do idsk=1,3                                                      3d25s20
        iaddk=nlk*(idsk-1+3*ispinks)                                    3d28s20
        do ispinbs=0,1                                                  3d25s20
         do idsb=1,3                                                    3d25s20
          iaddb=nlb*(idsb-1+3*ispinbs)                                  3d28s20
          if(ispinks.eq.ispinbs)then                                    3d25s20
c
c     recall matrices coming out of onep and derid have ket index first
c
           jbo=ip2eri(idsb,idsk)                                        3d25s20
           do ib=0,nlb-1                                                3d31s20
            do ik=0,nlk-1                                                3d25s20
             ikp=ik+iaddk                                                3d25s20
             ibp=iprim+iaddb+nlb6*ikp                                    3d25s20
             bc(ibp+ib)=bc(jbo+ik)                                      3d31s20
            end do                                                      3d25s20
            jbo=jbo+nlk                                                 3d31s20
           end do                                                       3d25s20
          end if                                                        3d25s20
         end do                                                         3d25s20
        end do                                                          3d25s20
       end do                                                           3d25s20
      end do                                                            3d25s20
c
c     iprim is now matrix in 4c ss basis.
c     now transform to 2c basis
c
      itmpr=ibcoff                                                      3d31s20
      itmpi=itmpr+nlb6*nlk2                                             3d31s20
      ibcoff=itmpi+nlb6*nlk2                                            3d31s20
      call enough('ssdvrc.  2',bc,ibc)
      do i=itmpr,ibcoff-1                                               3d31s20
       bc(i)=0d0                                                        3d31s20
      end do                                                            3d31s20
      do ispinkl=0,1                                                    3d31s20
       jspinkl=ispinkl+1                                                3d31s20
       do ik=0,nlk-1                                                    3d31s20
        ikp=ik+ispinkl*nlk                                              3d31s20
        iadr=itmpr+nlb6*ikp                                             3d31s20
        iadi=itmpi+nlb6*ikp                                             3d31s20
        do kk=1,6                                                       3d31s20
         iadf=iprim+nlb6*(ik+nlk*(kk-1))                                3d31s20
         do kb=0,nlb6-1                                                 3d31s20
          bc(iadr+kb)=bc(iadr+kb)+bc(iadf+kb)*cvec(kk,jspinkl,1)        3d31s20
          bc(iadi+kb)=bc(iadi+kb)+bc(iadf+kb)*cvec(kk,jspinkl,2)        3d31s20
         end do                                                         3d28s20
        end do                                                          3d31s20
       end do                                                           3d31s20
      end do                                                            3d31s20
      jtmpr=ibcoff                                                      3d31s20
      jtmpi=jtmpr+nlb6*nlk2                                             3d31s20
      ibcoff=jtmpi+nlb6*nlk2                                            3d31s20
      call enough('ssdvrc.  3',bc,ibc)
      do i=0,nlk2-1                                                     3d31s20
       do j=0,nlb6-1                                                    3d31s20
        ji=itmpr+j+nlb6*i                                               3d31s20
        ij=jtmpr+i+nlk2*j                                               3d31s20
        bc(ij)=bc(ji)                                                   3d31s20
        ji=itmpi+j+nlb6*i                                               3d31s20
        ij=jtmpi+i+nlk2*j                                               3d31s20
        bc(ij)=bc(ji)                                                   3d31s20
       end do                                                           3d31s20
      end do                                                            3d31s20
c
c     parah0 assumes I return the transpose
c
      do ispinbl=0,1                                                    3d31s20
       jspinbl=ispinbl+1                                                3d31s20
       do ib=0,nlb-1                                                    3d31s20
        ibp=ib+ispinbl*nlb                                              3d31s20
        iadr=iout+nlk2*ibp                                              3d31s20
        iadi=iouti+nlk2*ibp                                             3d31s20
        do kb=1,6                                                       3d31s20
         iadfr=jtmpr+nlk2*(ib+nlb*(kb-1))                               3d31s20
         iadfi=jtmpi+nlk2*(ib+nlb*(kb-1))                               3d31s20
         do kk=0,nlk2-1                                                 3d31s20
          bc(iadr+kk)=bc(iadr+kk)+bc(iadfr+kk)*cvec(kb,jspinbl,1)       3d31s20
     $                           +bc(iadfi+kk)*cvec(kb,jspinbl,2)       3d31s20
          bc(iadi+kk)=bc(iadi+kk)+bc(iadfi+kk)*cvec(kb,jspinbl,1)       3d31s20
     $                           -bc(iadfr+kk)*cvec(kb,jspinbl,2)       3d31s20
         end do                                                         3d31s20
        end do                                                          3d31s20
       end do                                                           3d31s20
      end do                                                            3d31s20
      ibcoff=iout                                                       3d4s20
      return                                                            3d4s20
      end                                                               3d4s20
