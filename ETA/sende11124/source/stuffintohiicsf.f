c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine stuffintohiicsf(g,isto,ncsft,ihcpy,vec,nps,ipointf,    7d12s19
     $     icall,iter,ntail,nvv,hiv,nlzzu,npsk,veclz,bc,ibc)            11d10s22
      implicit real*8 (a-h,o-z)
      logical lprint                                                    6d26s18
c
c     take hc and pick out pspace parts,
c     and stuff into decimated hamiltonian
c
      include "common.store"
      dimension g(ncsft,isto),vec(ncsft,isto),hiv(*),ipointf(*),        12d27s19
     $     veclz(*)                                                     12d27s19
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data ncall/0/                                                     9d5s19
      save ncall                                                        9d5s19
      ncall=ncall+1                                                     9d5s19
      lprint=.false.
      if(lprint)then                                                    7d12s19
       write(6,*)('Hi, my name is stuffintohiicsf '),isto,nps,nvv,icall,
     $     ncall
       write(6,*)('what we have for g: ')
       call prntm2(g,ncsft,isto,ncsft)
      end if                                                            7d12s19
      if(nlzzu.ne.0)then                                                12d27s19
       itmp=ibcoff                                                      12d27s19
       itmpr=itmp+(npsk+nvv)*isto                                       12d27s19
       ibcoff=itmpr+nps*isto                                            12d27s19
       call enough('stuffintohiicsf.  1',bc,ibc)
       jtmpr=itmpr-1-isto                                               12d27s19
       do i=1,isto                                                      12d27s19
        do j=1,nps                                                      12d27s19
         ii=ipointf(j)                                                  12d27s19
         iad=jtmpr+i+isto*j                                             12d27s19
         bc(iad)=g(ii,i)                                                12d27s19
        end do                                                          12d27s19
       end do                                                           12d27s19
       call dgemm('n','n',isto,npsk,nps,1d0,bc(itmpr),isto,veclz,nps,
     $      0d0,bc(itmp),isto,                                          12d27s19
     d' stuffintohiicsf.  1')
       ibcoff=itmpr                                                     12d27s19
       npsu=npsk                                                        12d27s19
      else                                                              12d27s19
       npsu=nps                                                         12d27s19
       itmp=ibcoff
       ibcoff=itmp+(nps+nvv)*isto                                        9d4s19
       call enough('stuffintohiicsf.  2',bc,ibc)
       do i=1,isto
        do j=1,nps
         ii=ipointf(j)                                                   7d12s19
         iad=itmp+i-1+isto*(j-1)
         bc(iad)=g(ii,i)                                                 7d12s19
        end do
       end do
      end if                                                            12d27s19
      if(nvv.gt.0)then                                                  9d4s19
       itmpp=itmp+isto*npsu                                             12d27s19
       igt=ibcoff                                                       9d4s19
       ibcoff=igt+isto*ncsft                                            9d4s19
       call enough('stuffintohiicsf.  3',bc,ibc)
       do i=1,isto                                                      9d4s19
        do j=1,ncsft                                                    9d4s19
         ij=igt+i-1+isto*(j-1)                                          9d4s19
         bc(ij)=vec(j,i)                                                9d5s19
        end do                                                          9d4s19
       end do                                                           9d4s19
       call dgemm('n','n',isto,nvv,ncsft,1d0,bc(igt),isto,hiv,ncsft,0d0,9d4s19
     $     bc(itmpp),isto,                                              9d4s19
     d' stuffintohiicsf.  2')
       ibcoff=igt                                                       9d4s19
      end if                                                            9d4s19
      if(lprint)then                                                    6d26s18
       write(6,*)('rectangle of p-space root hamiltonian matrix ')
       call prntm2(bc(itmp),isto,npsu+nvv,isto)                         12d27s19
      end if                                                            6d26s18
      ii=0
      nn=0
      do i=1,npsu+nvv                                                   12d27s19
       do j=1,i
        if(mod(ii,mynprocg).eq.mynowprog)then
         nn=nn+1
        end if
        ii=ii+1
       end do
      end do
      nnew=npsu+isto+nvv                                                12d27s19
      ntri=(nnew*(nnew+1))/2                                            6d27s18
      nhx=ntri/mynprocg                                                 6d27s18
      nleft=ntri-nhx*mynprocg                                           2d3s22
      if(mynowprog.lt.nleft)nhx=nhx+1                                   6d27s18
      ntail=nhx-nn+1                                                    6d27s18
      ihcpy=ibcoff                                                      6d27s18
      ibcoff=ihcpy+ntail                                                6d27s18
      call enough('stuffintohiicsf.  4',bc,ibc)
      xnan=-acos(-1d0)
      do i=0,ntail-1
       bc(ihcpy+i)=xnan
      end do
      jhcpy=ihcpy
      do i=1,isto
       ip=i+npsu+nvv                                                    12d27s19
       do j=1,npsu+nvv                                                  12d27s19
        iad=((ip*(ip-1))/2)+j-1
        if(mod(iad,mynprocg).eq.mynowprog)then                          6d22s18
         iad2=itmp+i-1+isto*(j-1)                                       7d12s19
         if(lprint)write(6,*)('storing '),ip,j,bc(iad2),iad                 6d26s18
         bc(jhcpy)=bc(iad2)                                             6d22s18
         jhcpy=jhcpy+1                                                  6d22s18
        end if                                                          6d22s18
       end do                                                           6d22s18
       do j=1,i                                                         6d22s18
        jp=j+npsu+nvv                                                   12d27s19
        iad=((ip*(ip-1))/2)+jp-1                                        6d22s18
        if(mod(iad,mynprocg).eq.mynowprog)then                          6d22s18
         dot=0d0                                                        6d22s18
         do l=1,ncsft                                                   6d22s18
          dot=dot+g(l,i)*vec(l,j)                                       7d12s19
         end do                                                         6d22s18
         bc(jhcpy)=dot                                                  6d22s18
         jhcpy=jhcpy+1                                                  6d22s18
        end if                                                          6d22s18
       end do                                                           6d22s18
      end do
      nigot=jhcpy-ihcpy
      return
      end
