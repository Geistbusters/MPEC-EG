c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine stuffintohii2(g,ncsft,vecqs,ihcpy,ncsfq,nrootw,savqv,  3d25s21
     $     savqg,savpg,ipointf,nps,ipointq,ntail,nvv,hiv,nlzzu,npsk,    3d25s21
     $     veclz,ngot,gq,nigot,iprint,bc,ibc)                           11d10s22
      implicit real*8 (a-h,o-z)
      logical lprint                                                    6d26s18
c
c     take hc and pick out pspace parts,
c     and stuff into decimated hamiltonian
c
      include "common.store"
      dimension g(ncsft,*),vecqs(ncsfq,*),hiv(*),ipointf(*),            3d25s21
     $     veclz(*),savqv(ncsfq,*),savqg(ncsfq,*),ipointq(*),           3d25s21
     $     savpg(nps,*),gq(ncsfq,*)                                     3d25s21
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data ncall/0/                                                     9d5s19
      save ncall                                                        9d5s19
      ncall=ncall+1                                                     9d5s19
      lprint=iprint.ne.0                                                3d25s21
      if(lprint)then                                                    7d12s19
       write(6,*)('Hi, my name is stuffintohii2 '),nrootw,nps,nvv,
     $     ncall,ngot
       write(6,*)('what we have for g: ')
       call prntm2(g,ncsft,nrootw,ncsft)
       write(6,*)('what we have for vecqs: ')
       call prntm2(vecqs,ncsfq,nrootw,ncsfq)                            3d25s21
       if(ngot.gt.0)then                                                3d25s21
        write(6,*)('what we have for savqv: ')
        call prntm2(savqv,ncsfq,ngot,ncsfq)
        write(6,*)('what we have for savqg: ')
        call prntm2(savqg,ncsfq,ngot,ncsfq)
        write(6,*)('what we have for savpg: ')
        call prntm2(savpg,nps,ngot,nps)
       end if                                                           3d25s21
      end if                                                            7d12s19
      isto=ngot+nrootw                                                  3d25s21
      if(nlzzu.ne.0)then                                                12d27s19
       itmp=ibcoff                                                      12d27s19
       itmpr=itmp+(npsk+nvv)*isto                                       3d25s21
       ibcoff=itmpr+nps*isto                                            3d25s21
       call enough('stuffintohii2.  1',bc,ibc)
       jtmpr=itmpr-1-isto                                               3d25s21
       do i=1,ngot                                                      3d25s21
        do j=1,nps                                                      3d25s21
         iad=jtmpr+i+isto*j                                             3d25s21
         bc(iad)=savpg(j,i)                                             3d25s21
        end do                                                          3d25s21
       end do                                                           3d25s21
       jtmpr=jtmpr+ngot                                                 3d25s21
       do i=1,nrootw                                                    3d25s21
        do j=1,nps                                                      12d27s19
         ii=ipointf(j)                                                  12d27s19
         iad=jtmpr+i+isto*j                                             3d25s21
         bc(iad)=g(ii,i)                                                12d27s19
        end do                                                          12d27s19
       end do                                                           12d27s19
       if(lprint)then                                                   3d25s21
        write(6,*)('ps part going into veclz ')
        call prntm2(bc(itmpr),isto,nps,isto)
       end if                                                           3d25s21
       if(isto.gt.0)then                                                3d28s22
        call dgemm('n','n',isto,npsk,nps,1d0,bc(itmpr),isto,veclz,       3d25s21
     $      nps,0d0,bc(itmp),isto,                                      3d25s21
     d' stuffintohii2.  1')
       end if                                                           3d28s22
       ibcoff=itmpr                                                     12d27s19
       npsu=npsk                                                        12d27s19
      else                                                              12d27s19
       npsu=nps                                                         12d27s19
       itmp=ibcoff
       ibcoff=itmp+(nps+nvv)*isto                                       3d25s21
       call enough('stuffintohii2.  2',bc,ibc)
       jtmp=itmp-1-isto                                                 3d25s21
       do i=1,ngot                                                      3d25s21
        do j=1,nps                                                      3d25s21
         iad=jtmp+i+isto*j                                              3d25s21
         bc(iad)=savpg(j,i)                                             3d25s21
        end do                                                          3d25s21
       end do                                                           3d25s21
       jtmp=jtmp+ngot                                                   3d25s21
       do i=1,nrootw                                                    3d25s21
        do j=1,nps
         ii=ipointf(j)                                                   7d12s19
         iad=jtmp+i+isto*j                                              3d25s21
         bc(iad)=g(ii,i)                                                 7d12s19
        end do
       end do
      end if                                                            12d27s19
      if(nvv.gt.0)then                                                  9d4s19
       itmpp=itmp+isto*npsu                                             3d25s21
       igt=ibcoff                                                       9d4s19
       ibcoff=igt+isto*ncsfq                                            3d25s21
       call enough('stuffintohii2.  3',bc,ibc)
       jgt=igt-1-isto                                                   3d25s21
       do i=1,ngot                                                      3d25s21
        do j=1,ncsfq                                                    3d25s21
         ij=jgt+i+isto*j                                                3d25s21
         bc(ij)=savqv(j,i)                                              3d25s21
        end do                                                          3d25s21
       end do                                                           3d25s21
       jgt=jgt+ngot                                                     3d25s21
       do i=1,nrootw                                                    3d25s21
        do j=1,ncsfq                                                    3d25s21
         ij=jgt+i+isto*j                                                3d25s21
         bc(ij)=vecqs(j,i)                                              3d25s21
        end do                                                          9d4s19
       end do                                                           9d4s19
       call dgemm('n','n',isto,nvv,ncsfq,1d0,bc(igt),isto,hiv,          3d25s21
     $      ncsfq,0d0,bc(itmpp),isto,                                   3d25s21
     d' stuffintohii2.  2')
       ibcoff=igt                                                       9d4s19
      end if                                                            9d4s19
      if(lprint)then                                                    6d26s18
       write(6,*)('rectangle of p-space root hamiltonian matrix ')
       call prntm2(bc(itmp),isto,npsu+nvv,isto)                         3d25s21
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
      nnew=npsu+isto+nvv                                                3d25s21
      ntri=(nnew*(nnew+1))/2                                            6d27s18
      nhx=ntri/mynprocg                                                 6d27s18
      nleft=ntri-nhx*mynprocg                                           1d19s23
      if(mynowprog.lt.nleft)nhx=nhx+1                                   6d27s18
      ntail=nhx-nn+1                                                    6d27s18
      ihcpy=ibcoff                                                      6d27s18
      ibcoff=ihcpy+ntail                                                6d27s18
      call enough('stuffintohii2.  4',bc,ibc)
      xnan=-acos(-1d0)
      do i=0,ntail-1
       bc(ihcpy+i)=xnan
      end do
      jhcpy=ihcpy
      do i=1,isto+ngot                                                  2d3s22
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
         dot1=0d0
         dot2=0d0
         dot3=0d0
         dot4=0d0
         if(i.le.ngot)then                                              3d25s21
          do l=1,ncsfq                                                  3d25s21
           dot=dot+savqg(l,i)*savqv(l,j)                                3d25s21
          end do                                                        3d25s21
          dot1=dot
         else                                                           3d25s21
          im=i-ngot                                                     3d25s21
          if(j.le.ngot)then                                             3d25s21
           do l=1,ncsfq                                                 3d25s21
            dot=dot+savqv(l,j)*gq(l,im)                                 3d25s21
           end do                                                       3d25s21
           dot2=dot
          else                                                          3d25s21
           jm=j-ngot                                                    3d25s21
           do l=1,ncsfq                                                 3d25s21
            dot=dot+vecqs(l,jm)*gq(l,im)                                3d25s21
            dot4=dot4+vecqs(l,jm)**2
           end do                                                       3d25s21
           dot3=dot
          end if                                                        3d25s21
         end if                                                         3d25s21
         bc(jhcpy)=dot                                                  6d22s18
         jhcpy=jhcpy+1                                                  6d22s18
        end if                                                          6d22s18
       end do                                                           6d22s18
      end do
      do ir=1,nrootw                                                    3d25s21
       irp=ir+ngot                                                      3d25s21
       do j=1,ncsfq                                                     3d25s21
        savqv(j,irp)=vecqs(j,ir)                                        3d25s21
        savqg(j,irp)=gq(j,ir)                                           3d25s21
       end do                                                           3d25s21
       do j=1,nps                                                       3d25s21
        savpg(j,irp)=g(ipointf(j),ir)                                   3d25s21
       end do                                                           3d25s21
      end do                                                            3d25s21
      if(lprint)then                                                    3d25s21
       write(6,*)('what is now under savqv')
       call prntm2(savqv,ncsfq,isto,ncsfq)
       call dgemm('t','n',isto,isto,ncsfq,1d0,savqv,ncsfq,savqv,ncsfq,
     $      0d0,bc(ibcoff),isto,
     d' stuffintohii2.  3')
       write(6,*)('ortho test ')
       call prntm2(bc(ibcoff),isto,isto,isto)
       write(6,*)('what is now under savqg')
       call prntm2(savqg,ncsfq,isto,ncsfq)
       write(6,*)('what is now under savpg')
       call prntm2(savpg,nps,isto,nps)
      end if                                                            3d25s21
      ngot=isto                                                         3d22s22
       nigot=jhcpy-ihcpy
      if(lprint)then                                                    3d25s21
       write(6,*)('my part of hcpy ')
       call prntm2(bc(ihcpy),1,nigot,1)
      end if                                                            3d25s21
      ibcoff=itmp
      return
      end
