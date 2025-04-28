c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine assqn(vecold,nbasisp,nbasisc,isb,nlzz,ihit,            12d8s22
     $             natom,ngaus,ibdat,isym,iapair,ibstor,                12d8s22
     $             isstor,idorel,ascale,multh,nsymb,iqnw1,iqnw2,bc,ibc, 9d19s23
     $             epssym)                                              9d19s23
      implicit real*8 (a-h,o-z)                                         12d7s22
      integer*8 ihit                                                    12d8s22
      dimension vecold(*),nbasisp(*),nbasisc(*),                        12d8s22
     $     multh(8,8),iorb(8),noc(8),ixmtr(8,6),iorbq(8),iorbqz(8),     12d8s22
     $     iqnw1(*),iqnw2(*),ihit(*),islz(3)                            2d3s23
      include "common.store"                                            12d7s22
      ibcoffo=ibcoff                                                    12d7s22
      ncomp=1                                                           12d4s22
      if(idorel.ne.0)ncomp=2                                            12d4s22
      ib4c=ibcoff                                                       12d4s22
      do jsb=1,nsymb                                                    12d4s22
       iorb(jsb)=ibcoff                                                 12d4s22
       ibcoff=iorb(jsb)+nbasisp(jsb)*ncomp*nbasisp(jsb)                 12d4s22
       noc(jsb)=0                                                       12d4s22
      end do                                                            12d4s22
      call enough('assqn.1',bc,ibc)
      do iz=ib4c,ibcoff-1                                               12d4s22
       bc(iz)=0d0                                                       12d4s22
      end do                                                            12d4s22
      do jsb=1,nsymb                                                    12d4s22
       iorbq(jsb)=ibcoff                                                12d4s22
       ibcoff=ibcoff+nbasisp(jsb)                                       12d4s22
       iorbqz(jsb)=ibcoff                                               12d4s22
       ibcoff=ibcoff+nbasisp(jsb)                                       12d4s22
      end do                                                            12d4s22
      nbasdws=nbasisp(isb)*ncomp                                        12d7s22
      noc(isb)=nbasisc(isb)
      jorb=iorb(isb)-1                                                  12d4s22
      call enough('assqn.2',bc,ibc)
      do i=1,nbasdws*nbasisc(isb)                                       12d7s22
       bc(jorb+i)=vecold(i)                                              12d7s22
      end do                                                            12d7s22
      call propcas(natom,ngaus,ibdat,isym,iapair,ibstor,isstor,idorel,  12d7s22
     $     ascale,multh,nbasisp,iorb,ixmtr,noc,islz,nlzz,0,iorbq,1,     1d18s23
     $      iorbqz,bc,ibc,epssym)                                       9d19s23
      itmpqn=ibcoff                                                     2d3s23
      ibcoff=itmpqn+noc(isb)                                            2d3s23
      call enough('assqn.tmpqn',bc,ibc)                                 2d3s23
      do i=0,noc(isb)-1
       if(nlzz.eq.2)then
        ibc(itmpqn+i)=ibc(iqnw1(isb)+i)                                 2d3s23
       else
        ibc(itmpqn+i)=ibc(iqnw1(isb)+i)+100*ibc(iqnw2(isb)+i)           2d3s23
       end if
      end do
      do i=1,noc(isb)                                                   12d8s22
       ihit(i)=0                                                        12d8s22
      end do                                                            12d8s22
      nhit=0                                                            2d3s23
      do i=0,noc(isb)-1                                                 12d8s22
       if(nlzz.eq.2)then                                                2d2s23
        igoal=ibc(iorbq(isb)+i)                                          12d8s22
       else                                                             2d2s23
        igoal=ibc(iorbqz(isb)+i)+100*ibc(iorbq(isb)+1)                  2d2s23
       end if                                                           2d2s23
       do j=0,noc(isb)-1                                                12d8s22
        jp=j+1                                                          12d8s22
        if(ihit(jp).eq.0)then                                           12d8s22
         if(ibc(itmpqn+j).eq.igoal)then                                 2d3s23
          ihit(jp)=i+1                                                  12d8s22
          nhit=nhit+1                                                   2d3s23
          go to 446                                                     12d8s22
         end if                                                         12d8s22
        end if                                                          12d8s22
       end do                                                           12d8s22
  446  continue                                                         12d8s22
      end do                                                            12d8s22
      do i=0,noc(isb)-1                                                 12d8s22
       ip=i+1                                                           12d8s22
      end do                                                            12d8s22
      if(nhit.ne.noc(isb))then                                          2d3s23
       write(6,*)('oops! in assqn, I was only able match '),nhit,       2d3s23
     $     (' out of '),noc(isb),(' functions!')                        2d3s23
       write(6,*)('   want got ')                                       2d3s23
       do i=0,noc(isb)-1                                                2d3s23
        if(nlzz.eq.2)then                                                2d2s23
         igoal=ibc(iorbq(isb)+i)                                          12d8s22
        else                                                             2d2s23
         igoal=ibc(iorbqz(isb)+i)+100*ibc(iorbq(isb)+1)                  2d2s23
        end if                                                           2d2s23
        write(6,*)i+1,ibc(itmpqn+i),igoal                               2d3s23
       end do                                                           2d3s23
       stop                                                             2d3s23
      end if                                                            2d3s23
      icpy=ibcoff                                                       12d8s22
      ibcoff=icpy+noc(isb)*nbasdws                                      12d8s22
      call enough('assqn.4',bc,ibc)
      jcpy=icpy-1                                                       12d8s22
      do i=0,noc(isb)-1                                                 12d8s22
       ip=i+1                                                           12d8s22
       iad=nbasdws*(ihit(ip)-1)                                         12d8s22
       do j=1,nbasdws                                                   12d8s22
        bc(jcpy+j)=vecold(j+iad)                                        12d8s22
       end do                                                           12d8s22
       jcpy=jcpy+nbasdws                                                12d8s22
      end do                                                            12d8s22
      jcpy=icpy-1                                                       12d8s22
      do j=1,noc(isb)*nbasdws                                           12d8s22
       vecold(j)=bc(jcpy+j)                                             12d8s22
      end do                                                            12d8s22
      ibcoff=ibcoffo
      return
      end
