c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine genr(term1,term2,lhomo,jspin,nspin,ism,aww,isinfo,     8d8s22
     $     nstate,maxst1,nec,eweight,maxst2,pspacex,pthrs,ispina)       6d16s23
      implicit real*8 (a-h,o-z)
      character*(*) term1,term2
      logical lhomo,lsame                                               7d30s22
      integer*8 jspin(*)                                                8d8s22
      parameter (id=6,ida=2,idx=300)
      character*1 atom(2),crr,cii,cri,cir                               10d23s15
      character*2 type(id)
      character*32 data,trial(idx),trialp(idx),datac
      character*6 sname                                                 8d8s22
      dimension ipow(3,5,id),nml(id),is(ida),e(ida),it(ida),ic(ida),
     $     ee(idx),isym(idx),iprog(2,idx),isort(idx),ndeg(idx),iss(idx),
     $     ml(5,id),iml(2,idx),mlp(5,id),irpl(4),impl(4),irpr(4),impr(4)
     $     ,ihit(idx),isinfo(11,*),eweight(maxst2,*),pspacex(*)         8d8s22
      nstat=0                                                           8d1s22
      type(1)='S '
      ipow(1,1,1)=0
      ipow(2,1,1)=0
      ipow(3,1,1)=0
      nml(1)=1
      ml(1,1)=0
      mlp(1,1)=1
      type(2)='So'
      ipow(1,1,2)=1
      ipow(2,1,2)=1
      ipow(3,1,2)=1
      nml(2)=1
      ml(1,2)=0
      mlp(1,2)=-1
      type(3)='P '
      ipow(1,1,3)=1
      ipow(2,1,3)=1
      ipow(3,1,3)=0
      ml(1,3)=0
      mlp(1,3)=-1
      ipow(1,2,3)=1
      ipow(2,2,3)=0
      ipow(3,2,3)=1
      ml(2,3)=1
      mlp(2,3)=3
      ipow(1,3,3)=0
      ipow(2,3,3)=1
      ipow(3,3,3)=1
      ml(3,3)=1
      mlp(3,3)=-2
      nml(3)=3
      type(4)='Po'
      ipow(1,1,4)=1
      ipow(2,1,4)=0
      ipow(3,1,4)=0
      ml(1,4)=1
      mlp(1,4)=2
      ipow(1,2,4)=0
      ipow(2,2,4)=1
      ipow(3,2,4)=0
      ml(2,4)=1
c                                                                       1d25s16
c     I don't remember what mlp represents, but pattern                 1d25s16
c     is 01 is always minus, except here                                1d25s16
c                                                                       1d25s16
      mlp(2,4)=-1                                                       1d25s16
      ipow(1,3,4)=0
      ipow(2,3,4)=0
      ipow(3,3,4)=1
      ml(3,4)=0
      mlp(3,4)=3
      nml(4)=3
      type(5)='D '
      ipow(1,1,5)=0
      ipow(2,1,5)=0
      ipow(3,1,5)=0
      ml(1,5)=0
      mlp(1,5)=1
      ipow(1,2,5)=2
      ipow(2,2,5)=2
      ipow(3,2,5)=2
      ml(2,5)=2
      mlp(2,5)=3
      ipow(1,3,5)=1
      ipow(2,3,5)=1
      ipow(3,3,5)=0
      ml(3,5)=2
      mlp(3,5)=-2
      ipow(1,4,5)=1
      ipow(2,4,5)=0
      ipow(3,4,5)=1
      ml(4,5)=1
      mlp(4,5)=5
      ipow(1,5,5)=0
      ipow(2,5,5)=1
      ipow(3,5,5)=1
      ml(5,5)=1
      mlp(5,5)=-4
      nml(5)=5
      type(6)='Do'
      ipow(1,1,6)=0
      ipow(2,1,6)=0
      ipow(3,1,6)=1
      ml(1,6)=2
      mlp(1,6)=4
      ipow(1,2,6)=1
      ipow(2,2,6)=0
      ipow(3,2,6)=0
      ml(2,6)=1
      mlp(2,6)=3
      ipow(1,3,6)=0
      ipow(2,3,6)=1
      ipow(3,3,6)=0
      ml(3,6)=1
      mlp(3,6)=-2
      ipow(1,4,6)=3
      ipow(2,4,6)=3
      ipow(3,4,6)=3
      ml(4,6)=2
      mlp(4,6)=-1
      ipow(1,5,6)=1
      ipow(2,5,6)=1
      ipow(3,5,6)=1
      ml(5,6)=0
      mlp(5,6)=-5
      nml(6)=5
      ntrial=0
      it(1)=0
      it(2)=0
      do i=1,id
       if(term1(1:2).eq.type(i))it(1)=i
       if(term2(1:2).eq.type(i))it(2)=i
      end do
      if(min(it(1),it(2)).eq.0)then
       write(6,*)('type not matched: ')
       if(it(1).eq.0)write(6,*)term1
       if(it(2).eq.0)write(6,*)term2
       stop 'genr'                                                      6d16s23
      end if
      if(lhomo)then                                                     7d30s22
       do il=1,nml(it(1))
        if(mlp(il,it(1)).gt.0)then
         if(ml(il,it(1)).gt.0)then
          ncoml=2
          if(mod(ml(il,it(1)),2).eq.0)then
           ipsl=2
          else
           ipsl=1
          end if
         else
          ncoml=1
         end if
         irpl(1)=ipow(1,il,it(1))
         irpl(2)=ipow(2,il,it(1))
         irpl(3)=ipow(3,il,it(1))
         irpl(4)=1
         if(mlp(il,it(1)).ne.il)then
          impl(1)=ipow(1,mlp(il,it(1)),it(1))
          impl(2)=ipow(2,mlp(il,it(1)),it(1))
          impl(3)=ipow(3,mlp(il,it(1)),it(1))
          impl(4)=1
         else
          impl(1)=0
          impl(2)=0
          impl(3)=0
          impl(4)=0
         end if
        else if(mlp(il,it(1)).eq.-il)then
         ncoml=1
         irpl(1)=0
         irpl(2)=0
         irpl(3)=0
         irpl(4)=0
         impl(1)=ipow(1,il,it(1))
         impl(2)=ipow(2,il,it(1))
         impl(3)=ipow(3,il,it(1))
         impl(4)=1
        else
         ncoml=0
        end if
        do jl=1,nml(it(2))
         if(mlp(jl,it(2)).gt.0)then
          if(ml(jl,it(2)).gt.0)then
           ncomr=2
           if(mod(ml(jl,it(2)),2).eq.0)then
            ipsr=2
           else
            ipsr=1
           end if
          else
           ncomr=1
          end if
          irpr(1)=ipow(1,jl,it(2))
          irpr(2)=ipow(2,jl,it(2))
          irpr(3)=ipow(3,jl,it(2))
          irpr(4)=1
          if(mlp(jl,it(2)).ne.jl)then
           impr(1)=ipow(1,mlp(jl,it(2)),it(2))
           impr(2)=ipow(2,mlp(jl,it(2)),it(2))
           impr(3)=ipow(3,mlp(jl,it(2)),it(2))
           impr(4)=1
          else
           impr(1)=0
           impr(2)=0
           impr(3)=0
           impr(4)=0
          end if
         else if(mlp(jl,it(2)).eq.-jl)then
          ncomr=1
          irpr(1)=0
          irpr(2)=0
          irpr(3)=0
          irpr(4)=0
          impr(1)=ipow(1,jl,it(2))
          impr(2)=ipow(2,jl,it(2))
          impr(3)=ipow(3,jl,it(2))
          impr(4)=1
         else
          ncomr=0
         end if
         do iq=1,ncoml
          if(ncoml.gt.1)then
           if(ipsl.eq.1)then
            fri=(-1d0)**(iq+1)
            fii=-1d0
            fmi=fri
           else
            fri=1d0
            fii=(-1d0)**iq
            fmi=-fii
           end if
          else
           fri=1d0
           fii=1d0
           fmi=fri                                                      1d19s23
          end if
          do jc=1,ncomr
           if(ncomr.gt.1)then
            if(ipsr.eq.1)then
             frj=(-1d0)**(jc+1)
             fij=-1d0
             fmj=frj
            else
             frj=1d0
             fij=(-1d0)**jc
             fmj=-fij
            end if
           else
            frj=1d0
            fij=1d0
            fmj=frj                                                     1d19s23
           end if
           msum=nint(fmi*ml(il,it(1))+fmj*ml(jl,it(2)))
           if(fri*frj.gt.0.1d0)then
            crr='+'
           else if(fri*frj.lt.-0.1d0)then
            crr='-'
           else
            crr=' '
           end if
           if(irpl(4)*irpr(4).eq.0)crr='s '
           if(fii*fij.gt.0.1d0)then
            cii='+'
           else if(fii*fij.lt.-0.1d0)then
            cii='-'
           else
            cii=' '
           end if
           if(impl(4)*impr(4).eq.0)cii='s '
           if(fri*fij.gt.0.1d0)then
            cri='+'
           else if(fri*fij.lt.-0.1d0)then
            cri='-'
           else
            cri=' '
           end if
           if(irpl(4)*impr(4).eq.0)cri='s '
           if(fii*frj.gt.0.1d0)then
            cir='+'
           else if(fii*frj.lt.-0.1d0)then
            cir='-'
           else
            cir=' '
           end if
           if(impl(4)*irpr(4).eq.0)cir='s '
           if(msum.ge.0)then                                            7d30s22
            write(data(31:32),64)msum                                   8d1s22
            nterm=0
            data(16:30)='               '
            if(crr.ne.'s'.and.cii.ne.'s')then
             nterm=2
             write(data(1:15),252)crr,(irpl(jq),irpr(jq),jq=1,3)
             write(data(16:30),252)cii,(impl(jq),impr(jq),jq=1,3)
            else if(crr.ne.'s')then
             nterm=1
             write(data(1:15),252)crr,(irpl(jq),irpr(jq),jq=1,3)
            else if(cii.ne.'s')then
             nterm=1
             write(data(1:15),252)cii,(impl(jq),impr(jq),jq=1,3)
            end if
            if(nterm.gt.0)then
             ipass=0
 2100        continue
              do jq=1,ntrial
               if(data.eq.trial(jq).or.data.eq.trialp(jq))then
                go to 33
               end if
              end do
              if(ipass.eq.0)then
               if(data(16:16).ne.' ')then
                datac(1:15)=data(16:30)
                datac(16:30)=data(1:15)
                datac(31:32)=data(31:32)
                data=datac
                ipass=1
                go to 2100
               end if
              end if
              ntrial=ntrial+1
              if(ntrial.gt.idx)stop 'idx:trial'
              trial(ntrial)=data
              if(data(1:1).eq.'+')then
               data(1:1)='-'
              else
               data(1:1)='+'
              end if
              if(data(16:16).eq.'+')then
               data(16:16)='-'
              else if(data(16:16).eq.'-')then
               data(16:16)='+'
              end if
              trialp(ntrial)=data
   33        continue
            end if
            nterm=0
            data(16:30)='               '
            if(cri.ne.'s'.and.cir.ne.'s')then
             nterm=2
             write(data(1:15),252)cri,(irpl(jq),impr(jq),jq=1,3)
             write(data(16:30),252)cir,(impl(jq),irpr(jq),jq=1,3)
            else if(cri.ne.'s')then
             nterm=1
             write(data(1:15),252)cri,(irpl(jq),impr(jq),jq=1,3)
            else if(cir.ne.'s')then
             nterm=1
             write(data(1:15),252)cir,(impl(jq),irpr(jq),jq=1,3)
            end if
            if(nterm.gt.0)then
             ipass=0
 2101        continue
              do jq=1,ntrial
               if(data.eq.trial(jq).or.data.eq.trialp(jq))then
                go to 35
               end if
              end do
              if(ipass.eq.0)then
               if(data(16:16).ne.' ')then
                datac(1:15)=data(16:30)
                datac(16:30)=data(1:15)
                datac(31:32)=data(31:32)                                7d30s22
                data=datac
                ipass=1
                go to 2101
               end if
              end if
              ntrial=ntrial+1
              if(ntrial.gt.idx)stop 'idx:trial'
              trial(ntrial)=data
              if(data(1:1).eq.'+')then
               data(1:1)='-'
              else
               data(1:1)='+'
              end if
              if(data(16:16).eq.'+')then
               data(16:16)='-'
              else if(data(16:16).eq.'-')then
               data(16:16)='+'
              end if
              trialp(ntrial)=data
   35        continue
            end if
c     msum ge 0
           end if
c     do jc
          end do
c     do iq
         end do
  252    format(4(a1,3(2i2,1x),5x),7i3)
c     do jl
        end do
c     do il
       end do
       if(ntrial.gt.0)then
        do jq=1,ntrial
         ihit(jq)=0
        end do
        nok=0
        do jq=1,ntrial
         if(ihit(jq).eq.0)then
          data=trial(jq)
          data(3:3)=trial(jq)(5:5)
          data(5:5)=trial(jq)(3:3)
          data(8:8)=trial(jq)(10:10)
          data(10:10)=trial(jq)(8:8)
          data(13:13)=trial(jq)(15:15)
          data(15:15)=trial(jq)(13:13)
          data(31:32)=trial(jq)(31:32)                                  7d30s22
          lsame=data(3:3).eq.data(5:5).and.
     $        data(8:8).eq.data(10:10).and.
     $        data(13:13).eq.data(15:15)
          if(data(16:16).ne.' ')then
           data(18:18)=trial(jq)(20:20)
           data(20:20)=trial(jq)(18:18)
           data(23:23)=trial(jq)(25:25)
           data(25:25)=trial(jq)(23:23)
           data(28:28)=trial(jq)(30:30)
           data(30:30)=trial(jq)(28:28)
           if(.not.lsame)then
            lsame=data(3:3).eq.data(20:20).and.
     $          data(5:5).eq.data(18:18).and.
     $          data(8:8).eq.data(25:25).and.
     $          data(10:10).eq.data(23:23).and.
     $          data(13:13).eq.data(30:30).and.
     $          data(15:15).eq.data(28:28)
           end if
          end if
          read(data(3:5),*)ixa,ixb
          ix=mod(ixa+ixb,2)
          read(data(6:10),*)iya,iyb
          iy=mod(iya+iyb,2)
          read(data(11:15),*)iza,izb
          iz=mod(iza+izb,2)
          if(term1.ne.term2)then
           lsame=.false.
           match=jq
          else
           do jk=1,ntrial
            if(data.eq.trial(jk).or.(
     $          data(1:15).eq.trial(jk)(16:30).and.
     $          data(16:30).eq.trial(jk)(1:15).and.
     $          data(31:32).eq.trial(jk)(31:32)))then                    7d30s22
             if(iz.eq.0)then
              crr='+'
             else
              crr='-'
             end if
             match=jk
             ihit(jk)=1
             go to 102
            end if
            if(data.eq.trialp(jk).or.
     $          (data(1:15).eq.trialp(jk)(16:30).and.
     $          data(16:30).eq.trialp(jk)(1:15).and.
     $          data(31:32).eq.trialp(jk)(31:32)))then                  7d30s22
             if(iz.eq.0)then
              crr='-'
             else
              crr='+'
             end if
             match=jk
             ihit(jk)=1
             go to 102
            end if
           end do
           write(6,*)('unable to find match!! ')
           write(6,*)data
           match=0
  102      continue
           if(crr.eq.'+')then
            iz0=0
           else
            iz0=1
           end if
           if(lsame)then                                                6d16s23
            itry=ispina
            if(mod(itry,2).ne.0)then                                     6d16s23
            else                                                         6d16s23
             iz0=1-iz0                                                   6d16s23
            end if                                                       6d16s23
           end if                                                        6d16s23
          end if                                                        6d16s23
          do iz=0,1
           if(lsame)then                                                7d30s22
            if(iz0.eq.iz)then                                           7d30s22
            else                                                        7d30s22
            end if                                                      7d30s22
           else                                                         7d30s22
            if(iz.eq.0)then                                             7d30s22
            else                                                        7d30s22
            end if
           end if                                                       7d30s22
           isum=ix+iy+iz
           if(isum.eq.0)then
            isymb=1
           else if(isum.eq.3)then
            isymb=8
           else if(isum.eq.1)then
            if(ix.eq.1)then
             isymb=2
            else if(iy.eq.1)then
             isymb=3
            else if(iz.eq.1)then
             isymb=5
            end if
           else
            if(ix.eq.0)then
             isymb=7
            else if(iy.eq.0)then
             isymb=6
            else
             isymb=4
            end if
           end if
  101      format(a32,' matches ',i3,3i2,' with phase ',a1,' lsame ',
     $        l1)
           if(lsame)then                                                8d1s22
            do js=1,nspin                                               8d1s22
             isum=ism+jspin(js)-3                                       8d1s22
             isum=isum/2                                                8d1s22
             if(mod(isum,2).eq.0.and.iz0.eq.iz)then                     8d9s22
              read(trial(jq)(31:32),*)lambda                            8d1s22
              do nn=1,nstate                                            8d8s22
               if(isinfo(1,nn).eq.nec.and.isinfo(2,nn).eq.jspin(js).and.8d8s22
     $            isinfo(3,nn).eq.isymb.and.isinfo(5,nn).eq.lambda)then 8d8s22
                isinfo(4,nn)=isinfo(4,nn)+1                             8d8s22
                if(isinfo(4,nn).gt.maxst2)then                          8d8s22
                 write(6,*)('too many roots!!! '),isinfo(4,nn),maxst2   8d8s22
                 write(6,*)('the parameter maxst2 in common.cas must'), 11d7s23
     $                (' be increased to at least '),isinfo(4,nn)       11d7s23
                 stop 'genr'                                            8d8s22
                end if                                                  8d8s22
                eweight(isinfo(4,nn),nn)=aww                            8d8s22
                go to 1492                                              8d1s22
               end if                                                   8d1s22
              end do                                                    8d1s22
              nstate=nstate+1                                             8d1s22
              if(nstate.gt.maxst1)then                                  8d8s22
               write(6,*)('too many states!!! '),nstate,maxst1          8d8s22
               write(6,*)('the parameter maxst1 in common.cas must'),   11d7s23
     $              (' be increased to at least '),nstate               11d7s23
               stop 'genr'                                              8d8s22
              end if                                                    8d8s22
              isinfo(1,nstate)=nec                                      8d8s22
              isinfo(2,nstate)=jspin(js)                                8d8s22
              isinfo(5,nstate)=lambda                                   8d8s22
              isinfo(3,nstate)=isymb                                    8d8s22
              isinfo(4,nstate)=1                                        8d8s22
              eweight(1,nstate)=aww                                     8d8s22
              pspacex(nstate)=pthrs                                     8d8s22
              if(lambda.eq.0)then                                       8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Sig'                                        8d8s22
               nx=nx+1                                                  8d8s22
               if(isymb.eq.1.or.isymb.eq.5)then                         8d8s22
                sname(nx:nx)='+'                                        8d8s22
               else                                                     8d8s22
                sname(nx:nx)='-'                                        8d8s22
               end if                                                   8d8s22
              else if(lambda.eq.1)then                                  8d8s22
               nx=2                                                     8d8s22
               sname(1:nx)='Pi'                                         8d8s22
              else if(lambda.eq.2)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Del'                                         8d8s22
              else if(lambda.eq.3)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Phi'                                         8d8s22
              else if(lambda.eq.4)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Gam'                                         8d8s22
              end if                                                    8d8s22
              nx=nx+1                                                   8d8s22
              if(isymb.eq.1.or.isymb.eq.4.or.                           8d8s22
     $           isymb.eq.6.or.isymb.eq.7)then                          8d8s22
               sname(nx:nx)='g'                                         8d8s22
              else                                                      8d8s22
               sname(nx:nx)='u'                                         8d8s22
              end if                                                    8d8s22
              do nn=1,nx                                                8d8s22
               isinfo(5+nn,nstate)=ichar(sname(nn:nn))                  8d8s22
              end do                                                    8d8s22
              do nn=nx+6,11                                             8d8s22
               isinfo(nn,nstate)=0                                      8d8s22
              end do                                                    8d8s22
 1492         continue
             else if(mod(isum,2).ne.0.and.iz0.ne.iz)then                8d9s22
              read(trial(jq)(31:32),*)lambda                             8d1s22
              do nn=1,nstate                                            8d8s22
               if(isinfo(1,nn).eq.nec.and.isinfo(2,nn).eq.jspin(js).and.8d8s22
     $            isinfo(3,nn).eq.isymb.and.isinfo(5,nn).eq.lambda)then 8d8s22
                isinfo(4,nn)=isinfo(4,nn)+1                             8d8s22
                if(isinfo(4,nn).gt.maxst2)then                          8d8s22
                 write(6,*)('too many roots!!! '),isinfo(4,nn),maxst2   8d8s22
                 write(6,*)('the parameter maxst2 in common.cas must'), 11d7s23
     $                (' be increased to at least '),isinfo(4,nn)       11d7s23
                 stop 'genr'                                            8d8s22
                end if                                                  8d8s22
                eweight(isinfo(4,nn),nn)=aww                            8d8s22
                go to 1493                                              8d1s22
               end if                                                   8d1s22
              end do                                                    8d1s22
              nstate=nstate+1                                             8d1s22
              if(nstate.gt.maxst1)then                                  8d8s22
               write(6,*)('too many states!!! '),nstate,maxst1          8d8s22
               write(6,*)('the parameter maxst1 in common.cas must'),   11d7s23
     $              (' be increased to at least '),nstate               11d7s23
               stop 'genr'                                              8d8s22
              end if                                                    8d8s22
              isinfo(1,nstate)=nec                                      8d8s22
              isinfo(2,nstate)=jspin(js)                                8d8s22
              isinfo(5,nstate)=lambda                                   8d8s22
              isinfo(3,nstate)=isymb                                    8d8s22
              isinfo(4,nstate)=1                                        8d8s22
              eweight(1,nstate)=aww                                     8d8s22
              pspacex(nstate)=pthrs                                     8d8s22
              if(lambda.eq.0)then                                       8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Sig'                                        8d8s22
               nx=nx+1                                                  8d8s22
               if(isymb.eq.1.or.isymb.eq.5)then                         8d8s22
                sname(nx:nx)='+'                                        8d8s22
               else                                                     8d8s22
                sname(nx:nx)='-'                                        8d8s22
               end if                                                   8d8s22
              else if(lambda.eq.1)then                                  8d8s22
               nx=2                                                     8d8s22
               sname(1:nx)='Pi'                                         8d8s22
              else if(lambda.eq.2)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Del'                                         8d8s22
              else if(lambda.eq.3)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Phi'                                         8d8s22
              else if(lambda.eq.4)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Gam'                                         8d8s22
              end if                                                    8d8s22
              nx=nx+1                                                   8d8s22
              if(isymb.eq.1.or.isymb.eq.4.or.                           8d8s22
     $           isymb.eq.6.or.isymb.eq.7)then                          8d8s22
               sname(nx:nx)='g'                                         8d8s22
              else                                                      8d8s22
               sname(nx:nx)='u'                                         8d8s22
              end if                                                    8d8s22
              do nn=1,nx                                                8d8s22
               isinfo(5+nn,nstate)=ichar(sname(nn:nn))                  8d8s22
              end do                                                    8d8s22
              do nn=nx+6,11                                             8d8s22
               isinfo(nn,nstate)=0                                      8d8s22
              end do                                                    8d8s22
 1493         continue
             end if                                                     8d1s22
            end do                                                      8d1s22
           else                                                         8d1s22
            read(trial(jq)(31:32),*)lambda                              8d1s22
            do js=1,nspin                                               8d1s22
              do nn=1,nstate                                            8d8s22
               if(isinfo(1,nn).eq.nec.and.isinfo(2,nn).eq.jspin(js).and.8d8s22
     $            isinfo(3,nn).eq.isymb.and.isinfo(5,nn).eq.lambda)then 8d8s22
                isinfo(4,nn)=isinfo(4,nn)+1                             8d8s22
                if(isinfo(4,nn).gt.maxst2)then                          8d8s22
                 write(6,*)('too many roots!!! '),isinfo(4,nn),maxst2   8d8s22
                 write(6,*)('the parameter maxst2 in common.cas must'), 11d7s23
     $                (' be increased to at least '),isinfo(4,nn)       11d7s23
                 stop 'genr'                                            8d8s22
                end if                                                  8d8s22
                eweight(isinfo(4,nn),nn)=aww                            8d8s22
               go to 1494                                               8d1s22
              end if                                                    8d1s22
             end do                                                     8d1s22
              nstate=nstate+1                                             8d1s22
              if(nstate.gt.maxst1)then                                  8d8s22
               write(6,*)('too many states!!! '),nstate,maxst1          8d8s22
               write(6,*)('the parameter maxst1 in common.cas must'),   11d7s23
     $              (' be increased to at least '),nstate               11d7s23
               stop 'genr'                                              8d8s22
              end if                                                    8d8s22
              isinfo(1,nstate)=nec                                      8d8s22
              isinfo(2,nstate)=jspin(js)                                8d8s22
              isinfo(5,nstate)=lambda                                   8d8s22
              isinfo(3,nstate)=isymb                                    8d8s22
              isinfo(4,nstate)=1                                        8d8s22
              eweight(1,nstate)=aww                                     8d8s22
              pspacex(nstate)=pthrs                                     8d8s22
              if(lambda.eq.0)then                                       8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Sig'                                        8d8s22
               nx=nx+1                                                  8d8s22
               if(isymb.eq.1.or.isymb.eq.5)then                         8d8s22
                sname(nx:nx)='+'                                        8d8s22
               else                                                     8d8s22
                sname(nx:nx)='-'                                        8d8s22
               end if                                                   8d8s22
              else if(lambda.eq.1)then                                  8d8s22
               nx=2                                                     8d8s22
               sname(1:nx)='Pi'                                         8d8s22
              else if(lambda.eq.2)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Del'                                         8d8s22
              else if(lambda.eq.3)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Phi'                                         8d8s22
              else if(lambda.eq.4)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Gam'                                         8d8s22
              end if                                                    8d8s22
              if(lhomo)then                                             8d9s22
               nx=nx+1                                                   8d8s22
               if(isymb.eq.1.or.isymb.eq.4.or.                           8d8s22
     $           isymb.eq.6.or.isymb.eq.7)then                          8d8s22
                sname(nx:nx)='g'                                         8d8s22
               else                                                      8d8s22
                sname(nx:nx)='u'                                         8d8s22
               end if                                                    8d8s22
              end if                                                    8d9s22
              do nn=1,nx                                                8d8s22
               isinfo(5+nn,nstate)=ichar(sname(nn:nn))                  8d8s22
              end do                                                    8d8s22
              do nn=nx+6,11                                             8d8s22
               isinfo(nn,nstate)=0                                      8d8s22
              end do                                                    8d8s22
 1494        continue                                                   8d1s22
            end do                                                      8d1s22
           end if                                                       8d1s22
          end do                                                        7d30s22
         end if
        end do
       end if
      else
       do il=1,nml(it(1))
        if(mlp(il,it(1)).gt.0)then
         if(ml(il,it(1)).gt.0)then
          ncoml=2
          if(mod(ml(il,it(1)),2).eq.0)then
           ipsl=2
          else
           ipsl=1
          end if
         else
          ncoml=1
         end if
         irpl(1)=ipow(1,il,it(1))
         irpl(2)=ipow(2,il,it(1))
         irpl(3)=ipow(3,il,it(1))
         irpl(4)=1
         if(mlp(il,it(1)).ne.il)then
          impl(1)=ipow(1,mlp(il,it(1)),it(1))
          impl(2)=ipow(2,mlp(il,it(1)),it(1))
          impl(3)=ipow(3,mlp(il,it(1)),it(1))
          impl(4)=1
         else
          impl(1)=0
          impl(2)=0
          impl(3)=0
          impl(4)=0
         end if
        else if(mlp(il,it(1)).eq.-il)then
         ncoml=1
         irpl(1)=0
         irpl(2)=0
         irpl(3)=0
         irpl(4)=0
         impl(1)=ipow(1,il,it(1))
         impl(2)=ipow(2,il,it(1))
         impl(3)=ipow(3,il,it(1))
         impl(4)=1
        else
         ncoml=0
        end if
        do jl=1,nml(it(2))
         if(mlp(jl,it(2)).gt.0)then
          if(ml(jl,it(2)).gt.0)then
           ncomr=2
           if(mod(ml(jl,it(2)),2).eq.0)then
            ipsr=2
           else
            ipsr=1
           end if
          else
           ncomr=1
          end if
          irpr(1)=ipow(1,jl,it(2))
          irpr(2)=ipow(2,jl,it(2))
          irpr(3)=ipow(3,jl,it(2))
          irpr(4)=1
          if(mlp(jl,it(2)).ne.jl)then
           impr(1)=ipow(1,mlp(jl,it(2)),it(2))
           impr(2)=ipow(2,mlp(jl,it(2)),it(2))
           impr(3)=ipow(3,mlp(jl,it(2)),it(2))
           impr(4)=1
          else
           impr(1)=0
           impr(2)=0
           impr(3)=0
           impr(4)=0
          end if
         else if(mlp(jl,it(2)).eq.-jl)then
          ncomr=1
          irpr(1)=0
          irpr(2)=0
          irpr(3)=0
          irpr(4)=0
          impr(1)=ipow(1,jl,it(2))
          impr(2)=ipow(2,jl,it(2))
          impr(3)=ipow(3,jl,it(2))
          impr(4)=1
         else
          ncomr=0
         end if
         do iq=1,ncoml
          if(ncoml.gt.1)then
           if(ipsl.eq.1)then
            fri=(-1d0)**(iq+1)
            fii=-1d0
            fmi=fri
           else
            fri=1d0
            fii=(-1d0)**iq
            fmi=-fii
           end if
          else
           fri=1d0
           fii=1d0
           fmi=fri                                                      1d19s23
          end if
          do jc=1,ncomr
           if(ncomr.gt.1)then
            if(ipsr.eq.1)then
             frj=(-1d0)**(jc+1)
             fij=-1d0
             fmj=frj
            else
             frj=1d0
             fij=(-1d0)**jc
             fmj=-fij
            end if
           else
            frj=1d0
            fij=1d0
            fmj=frj                                                     1d19s23
           end if
           msum=nint(fmi*ml(il,it(1))+fmj*ml(jl,it(2)))
           if(fri*frj.gt.0.1d0)then
            crr='+'
           else if(fri*frj.lt.-0.1d0)then
            crr='-'
           else
            crr=' '
           end if
           if(irpl(4)*irpr(4).eq.0)crr='s '
           if(fii*fij.gt.0.1d0)then
            cii='+'
           else if(fii*fij.lt.-0.1d0)then
            cii='-'
           else
            cii=' '
           end if
           if(impl(4)*impr(4).eq.0)cii='s '
           if(fri*fij.gt.0.1d0)then
            cri='+'
           else if(fri*fij.lt.-0.1d0)then
            cri='-'
           else
            cri=' '
           end if
           if(irpl(4)*impr(4).eq.0)cri='s '
           if(fii*frj.gt.0.1d0)then
            cir='+'
           else if(fii*frj.lt.-0.1d0)then
            cir='-'
           else
            cir=' '
           end if
           if(impl(4)*irpr(4).eq.0)cir='s '
           if(msum.ge.0)then                                            7d30s22
            write(data(31:32),64)msum                                    7d30s22
   64       format(i2)                                                   7d30s22
            nterm=0
            data(16:30)='               '
            if(crr.ne.'s'.and.cii.ne.'s')then
             nterm=2
             write(data(1:15),252)crr,(irpl(jq),irpr(jq),jq=1,3)
             write(data(16:30),252)cii,(impl(jq),impr(jq),jq=1,3)
            else if(crr.ne.'s')then
             nterm=1
             write(data(1:15),252)crr,(irpl(jq),irpr(jq),jq=1,3)
            else if(cii.ne.'s')then
             nterm=1
             write(data(1:15),252)cii,(impl(jq),impr(jq),jq=1,3)
            end if
            if(nterm.gt.0)then
             do jq=1,ntrial
              if(data.eq.trial(jq).or.data.eq.trialp(jq))then
               go to 2033
              end if
             end do
             ntrial=ntrial+1
             if(ntrial.gt.idx)stop 'idx:trial'
             trial(ntrial)=data
             if(data(1:1).eq.'+')then
              data(1:1)='-'
             else
              data(1:1)='+'
             end if
             if(data(16:16).eq.'+')then
              data(16:16)='-'
             else if(data(16:16).eq.'-')then
              data(16:16)='+'
             end if
             trialp(ntrial)=data
 2033        continue
            end if
            nterm=0
            data(16:30)='               '
            if(cri.ne.'s'.and.cir.ne.'s')then
             nterm=2
             write(data(1:15),252)cri,(irpl(jq),impr(jq),jq=1,3)
             write(data(16:30),252)cir,(impl(jq),irpr(jq),jq=1,3)
            else if(cri.ne.'s')then
             nterm=1
             write(data(1:15),252)cri,(irpl(jq),impr(jq),jq=1,3)
            else if(cir.ne.'s')then
             nterm=1
             write(data(1:15),252)cir,(impl(jq),irpr(jq),jq=1,3)
            end if
            if(nterm.gt.0)then
             do jq=1,ntrial
              if(data.eq.trial(jq).or.data.eq.trialp(jq))then
               go to 2035
              end if
             end do
             ntrial=ntrial+1
             if(ntrial.gt.idx)stop 'idx:trial'
             trial(ntrial)=data
             if(data(1:1).eq.'+')then
              data(1:1)='-'
             else
              data(1:1)='+'
             end if
             if(data(16:16).eq.'+')then
              data(16:16)='-'
             else if(data(16:16).eq.'-')then
              data(16:16)='+'
             end if
             trialp(ntrial)=data
 2035        continue
            end if
           end if                                                       7d30s22
          end do
         end do
        end do
       end do
       if(ntrial.gt.0)then
        do jq=1,ntrial
         read(trial(jq)(3:5),*)ixa,ixb
         ix=mod(ixa+ixb,2)
         read(trial(jq)(6:10),*)iya,iyb
         iy=mod(iya+iyb,2)
         read(trial(jq)(11:15),*)iza,izb
         iz=mod(iza+izb,2)
         if(ix.eq.1.and.iy.eq.1)then
          isymb=4
         else if(ix.eq.1)then
          isymb=2
         else if(iy.eq.1)then
          isymb=3
         else
          isymb=1
         end if
         read(trial(jq)(31:32),*)lambda
         do js=1,nspin                                                  8d1s22
              do nn=1,nstate                                            8d8s22
               if(isinfo(1,nn).eq.nec.and.isinfo(2,nn).eq.jspin(js).and.8d8s22
     $            isinfo(3,nn).eq.isymb.and.isinfo(5,nn).eq.lambda)then 8d8s22
                isinfo(4,nn)=isinfo(4,nn)+1                             8d8s22
                if(isinfo(4,nn).gt.maxst2)then                          8d8s22
                 write(6,*)('too many roots!!! '),isinfo(4,nn),maxst2   8d8s22
                 write(6,*)('the parameter maxst2 in common.cas must'), 11d7s23
     $                (' be increased to at least '),isinfo(4,nn)       11d7s23
                 stop 'genr'                                            8d8s22
                end if                                                  8d8s22
                eweight(isinfo(4,nn),nn)=aww                            8d8s22
           go to 1495                                                   8d1s22
          end if                                                        8d1s22
         end do                                                         8d1s22
              nstate=nstate+1                                             8d1s22
              if(nstate.gt.maxst1)then                                  8d8s22
               write(6,*)('too many states!!! '),nstate,maxst1          8d8s22
               write(6,*)('the parameter maxst1 in common.cas must'),   11d7s23
     $              (' be increased to at least '),nstate               11d7s23
               stop 'genr'                                              8d8s22
              end if                                                    8d8s22
              isinfo(1,nstate)=nec                                      8d8s22
              isinfo(2,nstate)=jspin(js)                                8d8s22
              isinfo(5,nstate)=lambda                                   8d8s22
              isinfo(3,nstate)=isymb                                    8d8s22
              isinfo(4,nstate)=1                                        8d8s22
              eweight(1,nstate)=aww                                     8d8s22
              pspacex(nstate)=pthrs                                     8d8s22
              if(lambda.eq.0)then                                       8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Sig'                                        8d8s22
               nx=nx+1                                                  8d8s22
               if(isymb.eq.1.or.isymb.eq.5)then                         8d8s22
                sname(nx:nx)='+'                                        8d8s22
               else                                                     8d8s22
                sname(nx:nx)='-'                                        8d8s22
               end if                                                   8d8s22
              else if(lambda.eq.1)then                                  8d8s22
               nx=2                                                     8d8s22
               sname(1:nx)='Pi'                                         8d8s22
              else if(lambda.eq.2)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Del'                                         8d8s22
              else if(lambda.eq.3)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Phi'                                         8d8s22
              else if(lambda.eq.4)then                                  8d8s22
               nx=3                                                     8d8s22
               sname(1:nx)='Gam'                                         8d8s22
              end if                                                    8d8s22
              if(lhomo)then
               nx=nx+1                                                   8d8s22
               if(isymb.eq.1.or.isymb.eq.4.or.                           8d8s22
     $            isymb.eq.6.or.isymb.eq.7)then                          8d8s22
                sname(nx:nx)='g'                                         8d8s22
               else                                                      8d8s22
                sname(nx:nx)='u'                                         8d8s22
               end if                                                    8d8s22
              end if
              do nn=1,nx                                                8d8s22
               isinfo(5+nn,nstate)=ichar(sname(nn:nn))                  8d8s22
              end do                                                    8d8s22
              do nn=nx+6,11                                             8d8s22
               isinfo(nn,nstate)=0                                      8d8s22
              end do                                                    8d8s22
 1495    continue                                                       8d1s22
         end do
        end do
       end if
      end if
      return
      end
