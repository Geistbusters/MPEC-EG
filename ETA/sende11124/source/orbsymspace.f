c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine orbsymspace(nsymb,iorbsym,idoub,iact,noc,nvirt,mlx,    4d5s21
     $     inbr,ipbr,bc,ibc)                                            11d9s22
      implicit real*8 (a-h,o-z)                                         4d5s21
      dimension iorbsym(*),idoub(*),iact(*),noc(*),nvirt(*)             4d5s21
      include "common.store"                                            4d5s21
      mlx=0                                                             4d5s21
      do isb=1,nsymb                                                    4d5s21
       if(min(idoub(isb),iact(isb)).gt.0)then                           4d5s21
        do i=0,idoub(isb)-1                                             4d5s21
         do j=0,iact(isb)-1                                             4d5s21
          jp=j+idoub(isb)                                               4d5s21
          if(ibc(iorbsym(isb)+i).eq.ibc(iorbsym(isb)+jp))then           4d5s21
           mlx=max(mlx,ibc(iorbsym(isb)+i))                             4d5s21
          end if                                                        4d5s21
         end do                                                         4d5s21
        end do                                                          4d5s21
       end if                                                           4d5s21
       if(min(nvirt(isb),noc(isb)).gt.0)then                            4d5s21
        do i=0,noc(isb)-1                                               4d5s21
         do j=0,nvirt(isb)-1                                            4d5s21
          jp=j+noc(isb)                                                 4d5s21
          if(ibc(iorbsym(isb)+i).eq.ibc(iorbsym(isb)+jp))then           4d5s21
           mlx=max(mlx,ibc(iorbsym(isb)+i))                             4d5s21
          end if                                                        4d5s21
         end do                                                         4d5s21
        end do                                                          4d5s21
       end if                                                           4d5s21
      end do                                                            4d5s21
      inbr=ibcoff                                                       4d5s21
      ipbr=inbr+mlx                                                     4d5s21
      ibcoff=ipbr+mlx                                                   4d5s21
      call enough('orbsymspace.  1',bc,ibc)
      do i=inbr,ibcoff-1                                                4d5s21
       ibc(i)=0                                                         4d5s21
      end do                                                            4d5s21
      nbrt=0                                                            4d5s21
      do isb=1,nsymb                                                    4d5s21
       do l=0,mlx-1                                                     4d5s21
        ibc(ipbr+l)=0                                                   4d5s21
       end do                                                           4d5s21
       if(min(idoub(isb),iact(isb)).gt.0)then                           4d5s21
        do i=0,idoub(isb)-1                                             4d5s21
         if(ibc(iorbsym(isb)+i).gt.0)then                               4d5s21
          jnbr=ipbr+ibc(iorbsym(isb)+i)-1                               4d5s21
          do j=0,iact(isb)-1                                             4d5s21
           jp=j+idoub(isb)                                               4d5s21
           if(ibc(iorbsym(isb)+i).eq.ibc(iorbsym(isb)+jp))then           4d5s21
            ibc(jnbr)=ibc(jnbr)+1                                       4d5s21
            nbrt=nbrt+1                                                 4d5s21
           end if                                                        4d5s21
          end do                                                         4d5s21
         end if                                                         4d5s21
        end do                                                          4d5s21
       end if                                                           4d5s21
       if(min(nvirt(isb),noc(isb)).gt.0)then                            4d5s21
        do i=0,noc(isb)-1                                               4d5s21
         if(ibc(iorbsym(isb)+i).gt.0)then                               4d5s21
          jnbr=ipbr+ibc(iorbsym(isb)+i)-1                               4d5s21
          do j=0,nvirt(isb)-1                                            4d5s21
           jp=j+noc(isb)                                                 4d5s21
           if(ibc(iorbsym(isb)+i).eq.ibc(iorbsym(isb)+jp))then           4d5s21
            ibc(jnbr)=ibc(jnbr)+1                                       4d5s21
            nbrt=nbrt+1                                                 4d5s21
           end if                                                        4d5s21
          end do                                                         4d5s21
         end if                                                         4d5s21
        end do                                                          4d5s21
       end if                                                           4d5s21
       do l=0,mlx-1                                                     4d5s21
        if(ibc(ipbr+l).gt.ibc(inbr+l))then                              4d5s21
         ibc(inbr+l)=ibc(ipbr+l)                                        4d5s21
        end if                                                          4d5s21
       end do                                                           4d5s21
      end do                                                            4d5s21
      do l=0,mlx-1                                                      4d5s21
       ibc(ipbr+l)=ibcoff                                               4d5s21
       ibcoff=ibcoff+ibc(inbr+l)                                        4d5s21
      end do                                                            4d5s21
      call enough('orbsymspace.  2',bc,ibc)
      return                                                            4d5s21
      end                                                               4d5s21
