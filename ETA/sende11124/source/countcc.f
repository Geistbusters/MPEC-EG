c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine countcc(nsymo,nfdat,nsymb,nvirt,mdoub,ndoub,multh)     6d18s22
      implicit real*8 (a-h,o-z)                                         6d18s22
      dimension nfdat(5,4,*),nvirt(*),multh(8,8)                        6d18s22
      mdoub=0                                                           6d18s22
      ndoub=0                                                           6d18s22
      do isb=1,nsymb                                                    6d18s22
       isbv12=multh(isb,nsymo)                                          6d18s22
       neq=0                                                            6d18s22
       nnot=0                                                           6d18s22
       do isbv1=1,nsymb                                                 6d18s22
        isbv2=multh(isbv1,isbv12)                                       6d18s22
        if(isbv2.ge.isbv1)then                                          6d18s22
         if(isbv1.eq.isbv2)then                                         6d18s22
          neq=neq+nvirt(isbv1)                                          6d18s22
          nvv=(nvirt(isbv1)*(nvirt(isbv1)-1))/2                         6d18s22
         else                                                           6d18s22
          nvv=nvirt(isbv1)*nvirt(isbv2)                                 6d18s22
         end if                                                         6d18s22
         nnot=nnot+nvv                                                  6d18s22
        end if                                                          6d18s22
       end do                                                           6d18s22
       nsp=neq+nnot                                                     6d18s22
       mdoub=mdoub+nsp*nfdat(2,1,isb)                                   6d18s22
       ndoub=ndoub+nsp*nfdat(3,1,isb)                                   6d18s22
       do l=2,4                                                         6d18s22
        mdoub=mdoub+nnot*nfdat(2,l,isb)                                 6d18s22
        ndoub=ndoub+nnot*nfdat(3,l,isb)                                 6d18s22
       end do                                                           6d18s22
      end do                                                            6d18s22
      return                                                            6d18s22
      end                                                               6d18s22
