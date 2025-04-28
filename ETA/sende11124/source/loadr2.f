c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine loadr2(ibcode,iorb,nsymb,nbasdws,iorbao,ibcodex,bc,ibc)11d9s22
      implicit real*8 (a-h,o-z)
      dimension ibcode(*),iorb(*),nbasdws(*),iorbao(*)                  2d8s20
      logical ldebug
      include "common.basis"
      include "common.input"                                            11d20s19
      include "common.store"                                            11d20s19
      include "common.print"                                            4d15s24
      ldebug=iprtr(32).ne.0                                             4d15s24
      if(ldebug)write(6,*)('in loadr2 '),nsymb                          10d27s20
      do isb=1,nsymb                                                    1d10s19
       write(6,*)isb,nbasdws(isb)                                       2d8s20
       ibcode(isb)=ibcoff                                               1d10s19
       ibcoff=ibcode(isb)+nbasdws(isb)                                  2d8s20
       call enough('loadr2.  1',bc,ibc)
       read(1)(ibc(ibcode(isb)+i),i=0,nbasdws(isb)-1)                   2d8s20
       if(ldebug)                                                       10d27s20
     $    write(6,*)('ibcode: '),(ibc(ibcode(isb)+i),i=0,nbasdws(isb)-1)10d27s20
       do i=0,nbasdws(isb)-1                                            2d8s20
        if(ibc(ibcode(isb)+i).gt.ibcodex)ibcodex=ibc(ibcode(isb)+i)     2d8s20
       end do                                                           2d8s20
      end do                                                            1d10s19
      if(idorel.eq.0)then                                               2d8s20
       ncomp=1                                                          2d8s20
      else                                                              2d8s20
       ncomp=2                                                          2d8s20
      end if                                                            2d8s20
      do isb=1,nsymb                                                    1d10s19
       iorbao(isb)=ibcoff                                               2d8s20
       ibcoff=iorbao(isb)+nbasb(isb)*nbasdws(isb)*ncomp                 2d8s20
       call enough('loadr2.  2',bc,ibc)
       read(1)(bc(iorbao(isb)+i),i=0,nbasb(isb)*nbasdws(isb)*ncomp-1)   2d8s20
       iorb(isb)=ibcoff                                                 2d8s20
       ibcoff=iorb(isb)+nbasdws(isb)*nbasc(isb)                         2d9s20
       read(1)(bc(iorb(isb)+i),i=0,nbasdws(isb)*nbasc(isb)-1)           2d9s20
       if(ldebug)then
       write(6,*)('ob vectors for symmetry block '),isb                 1d10s19
       call prntm2(bc(iorb(isb)),nbasc(isb),nbasdws(isb),nbasc(isb))    2d9s20
       call dgemm('t','n',nbasdws(isb),nbasdws(isb),nbasc(isb),1d0,     2d9s20
     $      bc(iorb(isb)),nbasc(isb),bc(iorb(isb)),nbasc(isb),0d0,      2d9s20
     $      bc(ibcoff),nbasdws(isb),                                    2d8s20
     d' loadr2.  1')
       write(6,*)('ortho check ')
       call prntm2(bc(ibcoff),nbasdws(isb),nbasdws(isb),nbasdws(isb))   2d8s20
       end if                                                           10d27s20
      end do                                                            1d10s19
      return
      end
