c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine epdvec(nroot,ndm,vec,vecx,ncsft,ipointf,basisv,nx,nps, 9d4s19
     $     nvv,nlzzu,npsk,veclz,bc,ibc)                                 11d14s22
      implicit real*8 (a-h,o-z)                                         7d11s19
      dimension vec(ndm,nroot),vecx(ncsft,nroot),ipointf(*),            7d11s19
     $     basisv(ncsft,*),veclz(*)                                     12d27s19
      include "common.store"
c
c     from ps vectors, generate best overall vector
c
      if(nlzzu.ne.0)then                                                12d27s19
       itmp=ibcoff                                                      12d27s19
       ibcoff=itmp+nps*nroot                                            12d27s19
       call enough('epdvec.  1',bc,ibc)
       call dgemm('n','n',nps,nroot,npsk,1d0,veclz,nps,vec,ndm,0d0,     12d27s19
     $      bc(itmp),nps,                                               12d27s19
     d' epdvec.  1')
       do i=1,nroot                                                      7d11s19
        do j=1,ncsft                                                     7d11s19
         vecx(j,i)=0d0                                                   7d11s19
        end do                                                           7d11s19
        jtmp=itmp-1+nps*(i-1)                                           12d27s19
        do j=1,nps                                                       7d11s19
         vecx(ipointf(j),i)=bc(jtmp+j)                                  12d27s19
        end do                                                           7d11s19
       end do                                                            7d11s19
       ibcoff=itmp                                                      12d27s19
      else                                                              12d27s19
       do i=1,nroot                                                      7d11s19
        do j=1,ncsft                                                     7d11s19
         vecx(j,i)=0d0                                                   7d11s19
        end do                                                           7d11s19
        do j=1,nps                                                       7d11s19
         vecx(ipointf(j),i)=vec(j,i)                                     7d11s19
        end do                                                           7d11s19
       end do                                                            7d11s19
      end if                                                            12d27s19
      if(nx.gt.0)then                                                   7d12s19
       if(nlzzu.ne.0)then                                               12d27s19
        npsp=npsk+1+nvv                                                 12d27s19
       else                                                             12d27s19
        npsp=nps+1+nvv                                                   9d4s19
       end if                                                           12d27s19
       call dgemm('n','n',ncsft,nroot,nx,1d0,basisv,ncsft,vec(npsp,1),  7d12s19
     $      ndm,1d0,vecx,ncsft,
     d' epdvec.  2')
      end if                                                            7d12s19
      return                                                            7d11s19
      end                                                               7d11s19
