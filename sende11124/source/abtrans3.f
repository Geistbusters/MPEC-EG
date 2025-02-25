      subroutine abtrans3(dd3x,dd3xin,is1,is2,is3,is4,mcol2,noc,        4d15s24
     $     nbasisp,nvirt,ncomp,iorb,nrow,bc,ibc)                        4d8s24
      implicit real*8 (a-h,o-z)                                         4d8s24
      include "common.store"                                            4d8s24
      dimension dd3x(*),dd3xin(*),noc(*),nbasisp(*),nvirt(*),iorb(*)    4d8s24
      nbp2=nbasisp(is2)*ncomp
      nbp1=nbasisp(is1)*ncomp
      itmp2=ibcoff                                                      4d8s24
      itmp3=itmp2+nvirt(is1)*nbp2*mcol2                                 4d8s24
      ibcoff=itmp3+nbp2*nvirt(is2)
      call enough('dorb4v.t2',bc,ibc)
      if(is1.eq.is2)then                                                4d8s24
       do i34=0,mcol2-1                                                 4d8s24
        do i2=0,nvirt(is2)-1                                            4d8s24
         do i1=0,i2                                                     4d8s24
          itri=((i2*(i2+1))/2)+i1
          iad1=1+i34+mcol2*(i1+nvirt(is1)*i2)
          iad2=1+itri+nrow*i34
          dd3x(iad1)=dd3xin(iad2)
         end do
         do i1=i2+1,nvirt(is1)-1
          iad1=1+i34+mcol2*(i1+nvirt(is1)*i2)
          dd3x(iad1)=0d0
         end do
        end do
       end do
      else                                                              4d8s24
       do i34=0,mcol2-1                                                 4d8s24
        do i2=0,nvirt(is2)-1                                            4d8s24
         do i1=0,nvirt(is1)-1                                           4d8s24
          iad1=1+i34+mcol2*(i1+nvirt(is1)*i2)
          iad2=1+i1+nvirt(is1)*(i2+nvirt(is2)*i34)
          dd3x(iad1)=dd3xin(iad2)
         end do
        end do
       end do
      end if                                                            4d8s24
      do iv=0,nvirt(is2)-1                                              4d8s24
       ivp=iv+noc(is2)                                                  4d8s24
       do id=0,nbp2-1                                                   4d8s24
        idv=iorb(is2)+id+nbp2*ivp                                       4d8s24
        ivd=itmp3+iv+nvirt(is2)*id                                      4d8s24
        bc(ivd)=bc(idv)                                                 4d8s24
       end do                                                           4d8s24
      end do                                                            4d8s24
      nnn=nvirt(is1)*mcol2
      call dgemm('n','n',nnn,nbp2,nvirt(is2),1d0,
     $         dd3x,nnn,bc(itmp3),nvirt(is2),0d0,
     $         bc(itmp2),nnn,'dorb4v.2')
      mcol3=mcol2*nbasisp(is2)
      itmp4=itmp3+mcol3*nvirt(is1)
      itmp5=itmp4+nvirt(is1)*nbasisp(is1)
      ibcoff=itmp5+mcol3*nbasisp(is1)
      call enough('dorb4v.35',bc,ibc)
      do in=0,ncomp-1
       do i2=0,nbasisp(is2)-1
        i2p=i2+in*nbasisp(is2)
        do i1=0,nvirt(is1)-1
         i12=itmp2+mcol2*(i1+nvirt(is1)*i2p)
         i21=itmp3+mcol2*(i2+nbasisp(is2)*i1)
         do i34=0,mcol2-1
          bc(i21+i34)=bc(i12+i34)
         end do
        end do
       end do
       do iv=0,nvirt(is1)-1
        ivp=iv+noc(is1)
        do id=0,nbasisp(is1)-1
         idp=id+in*nbasisp(is1)
         iad1=itmp4+iv+nvirt(is1)*id
         iad2=iorb(is1)+idp+nbp1*ivp
         bc(iad1)=bc(iad2)
        end do
       end do
       nnn=mcol2*nbasisp(is2)
       call dgemm('n','n',nnn,nbasisp(is1),nvirt(is1),1d0,
     $          bc(itmp3),nnn,bc(itmp4),nvirt(is1),0d0,
     $          bc(itmp5),nnn,'dorb4x.1')
       if(is1.ne.is2)then                                               4d9s24
        do i34=0,mcol2-1
         do i2=0,nbasisp(is2)-1
          do i1=0,nbasisp(is1)-1
           iad1=itmp5+i34+mcol2*(i2+nbasisp(is2)*i1)
           iad2=1+i1+nbasisp(is1)*(i2+nbasisp(is2)*(i34
     $            +mcol2*in))
           dd3x(iad2)=bc(iad1)
          end do
         end do
        end do
       else                                                             4d9s24
        nrown=(nbasisp(is1)*(nbasisp(is1)+1))/2                         4d9s24
        do i34=0,mcol2-1
         do i2=0,nbasisp(is2)-1
          itri=1+((i2*(i2+1))/2)+nrown*(i34+mcol2*in)
          do i1=0,i2-1
           iad21=itmp5+i34+mcol2*(i2+nbasisp(is2)*i1)
           iad12=itmp5+i34+mcol2*(i1+nbasisp(is2)*i2)
           dd3x(itri+i1)=bc(iad21)+bc(iad12)                            4d9s24
          end do                                                        4d9s24
          iad22=itmp5+i34+mcol2*i2*(nbasisp(is2)+1)                     4d9s24
          dd3x(itri+i2)=bc(iad22)
         end do
        end do
       end if                                                           4d9s24
      end do
      ibcoff=itmp2
      return
      end
