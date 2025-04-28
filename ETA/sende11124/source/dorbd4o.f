c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dorbd4o(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,      6d24s24
     $     nsdlk1,isblk1,id4o,ioooo,ionex,bc,ibc,itt,nbasdws,igoal)           6d14s24
      implicit real*8 (a-h,o-z)
      dimension idoubo(*),irefo(*),noc(*),nvirt(*),                     6d24s24
     $     isblk(4,*),isblk1(4,*),id4o(*),ioooo(*),ionex(*),            6d14s24
     $     itt(*),nbasdws(*)                                            10d19s23
      include "common.store"                                            7d21s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      ihit=ibcoff                                                       7d28s23
      ibcoff=ihit+nsdlk                                                 7d28s23
      call enough('dorbd4o.ihit',bc,ibc)
      jhit=ihit-1                                                       7d28s23
      do is=1,nsdlk                                                     7d28s23
       ibc(jhit+is)=0                                                   7d28s23
      end do                                                            7d28s23
      do is=1,nsdlk                                                     7d28s23
       if(isblk(1,is).eq.isblk(2,is))then                               8d6s24
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              8d6s24
        iswitch=0                                                       8d6s24
       else                                                             8d6s24
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      8d6s24
        iswitch=1                                                       8d6s24
       end if                                                           8d6s24
       if(isblk(3,is).eq.isblk(4,is))then                               8d6s24
        ncol=(irefo(isblk(3,is))*(irefo(isblk(3,is))+1))/2              8d6s24
        jswitch=0                                                       8d6s24
       else                                                             8d6s24
        ncol=irefo(isblk(3,is))*irefo(isblk(4,is))                      8d6s24
        jswitch=1                                                       8d6s24
       end if                                                           8d6s24
       if(ibc(jhit+is).eq.0.and.min(nrow,ncol).gt.0)then                8d6s24
        i43=0                                                            7d28s23
        i21=0                                                            7d28s23
        i2143=0                                                          7d28s23
        i31=0                                                            7d28s23
        i3132=0                                                          7d28s23
        i3132b=0                                                         7d28s23
        i3132c=0                                                         7d28s23
        do js=is,nsdlk                                                  7d28s23
         if(isblk(1,is).eq.isblk(1,js).and.isblk(2,is).eq.isblk(2,js)    7d28s23
     $       .and.isblk(3,is).eq.isblk(4,js))i43=js                     7d28s23
         if(isblk(1,is).eq.isblk(3,js).and.isblk(2,is).eq.isblk(4,js)    7d28s23
     $       .and.isblk(3,is).eq.isblk(1,js))i31=js                     7d28s23
         if(isblk(1,is).eq.isblk(4,js).and.isblk(2,is).eq.isblk(3,js)    7d28s23
     $       .and.isblk(3,is).eq.isblk(1,js))i3132c=js                  7d28s23
        end do                                                           7d28s23
        if(min(i43,i31,i3132c).eq.0)then
         write(6,*)('for is '),is,i43,i31,i3132c                         7d28s23
         write(6,*)('we can not match perms! '),is,
     $        (isblk(j,is),j=1,4)
         call dws_synca
         call dws_finalize
         stop
        end if
        if(min(iswitch,jswitch).eq.1.and.min(ncol,nrow).gt.0)then       7d28s23
         do i4=0,irefo(isblk(4,is))-1                                   7d28s23
          do i3=0,irefo(isblk(3,is))-1                                   7d28s23
           do i2=0,irefo(isblk(2,is))-1                                   7d28s23
            do i1=0,irefo(isblk(1,is))-1                                   7d28s23
             iad1=id4o(is)+i1+irefo(isblk(1,is))*(i2+irefo(isblk(2,is))  7d28s23
     $            *(i3+irefo(isblk(3,is))*i4))                          7d28s23
             iad4=id4o(i43)+i1+irefo(isblk(1,is))*(i2                   7d28s23
     $            +irefo(isblk(2,is))*(i4+irefo(isblk(4,is))*i3))       7d28s23
             iad5=id4o(i31)+i3+irefo(isblk(3,is))*(i4+irefo(isblk(4,is))7d28s23
     $            *(i1+irefo(isblk(1,is))*i2))                          7d28s23
             iad8=id4o(i3132c)+i3+irefo(isblk(3,is))*(i4                7d28s23
     $            +irefo(isblk(4,is))*(i2+irefo(isblk(2,is))*i1))       7d28s23
             avg=(bc(iad1)+bc(iad4)+bc(iad5)+bc(iad8))*0.25d0           7d28s23
             bc(iad1)=avg                                               7d28s23
             bc(iad4)=avg                                               7d28s23
             bc(iad5)=avg                                               7d28s23
             bc(iad8)=avg                                               7d28s23
            end do                                                      7d28s23
           end do                                                       7d28s23
          end do                                                        7d28s23
         end do                                                         7d28s23
        end if                                                          7d28s23
        ibc(jhit+is)=1                                                  7d28s23
        ibc(jhit+i43)=1                                                 7d28s23
        ibc(jhit+i31)=1                                                 7d28s23
        ibc(jhit+i3132c)=1                                              7d28s23
       end if                                                           7d28s23
      end do                                                            7d28s23
      ibcoff=ihit                                                       7d28s23
      igoul=itt(1)+1+7*6                                                10d27s23
      do is=1,nsdlk                                                     10d27s23
       if(min(irefo(isblk(1,is)),irefo(isblk(2,is)),irefo(isblk(3,is)), 10d27s23
     $      irefo(isblk(4,is))).gt.0)then                               10d27s23
        isb1=isblk(1,is)                                                10d27s23
        isb2=isblk(2,is)                                                10d27s23
        isb=isblk(3,is)                                                 10d27s23
        jsb=isblk(4,is)                                                 10d27s23
        do is1=1,nsdlk1                                                  10d27s23
         if(isblk(1,is).eq.isblk1(1,is1)                                10d27s23
     $        .and.isblk(2,is).eq.isblk1(2,is1)                         10d27s23
     $        .and.isblk(3,is).eq.isblk1(3,is1))then                    10d27s23
          call ilimts(noc(isb),noc(jsb),mynprocg,mynowprog,il,ih,i1s,   10d27s23
     $        i1e,i2s,i2e)                                              10d27s23
          nhere=ih+1-il                                                   10d27s23
          call ilimts(noc(isb),nvirt(jsb),mynprocg,mynowprog,jl,jh,j1s,  10d27s23
     $        j1e,j2s,j2e)                                                  10d27s23
          mhere=jh+1-jl                                                   10d27s23
          if(max(nhere,mhere).gt.0)then                                 10d27s23
           if(isb.eq.jsb)then                                              10d27s23
            nrow=(noc(isb1)*(noc(isb1)+1))/2                             10d27s23
            nrowd=(irefo(isb1)*(irefo(isb1)+1))/2                        10d27s23
            ncol=(irefo(isb)*(irefo(isb)+1))/2                          10d27s23
            isw=0                                                          10d27s23
           else                                                            10d27s23
            nrow=noc(isb1)*noc(isb2)                                     10d27s23
            nrowd=irefo(isb1)*irefo(isb2)                                10d27s23
            ncol=irefo(isb)*irefo(jsb)                                  10d27s23
            isw=1                                                          10d27s23
           end if                                                         10d27s23
           nn=irefo(jsb)*irefo(isb)                                     10d27s23
           nmul=nrowd*irefo(isb)                                         10d27s23
           itmpd=ibcoff                                                  10d27s23
           ibcoff=itmpd+nmul*irefo(jsb)                                  10d27s23
           call enough('dorbd4o.tmpd0',bc,ibc)                           10d27s23
           do i4=0,irefo(jsb)-1                                          10d27s23
            if(isw.eq.0)then                                             10d27s23
             iadd=id4o(is)+nrowd*((i4*(i4+1))/2)                         10d27s23
             do i3=0,i4                                                  10d27s23
              ff=1d0                                                     10d27s23
              if(i3.ne.i4)ff=0.5d0                                       10d27s23
              jadd=iadd+nrowd*i3                                        10d27s23
              jtmpd1=itmpd+i3+irefo(jsb)*i4                             10d27s23
              jtmpd2=itmpd+i4+irefo(jsb)*i3                             10d27s23
              do i12=0,nrowd-1                                           10d27s23
               bc(jtmpd1)=bc(jadd+i12)*ff                                10d27s23
               bc(jtmpd2)=bc(jadd+i12)*ff                                10d27s23
               jtmpd1=jtmpd1+nn                                         10d27s23
               jtmpd2=jtmpd2+nn                                         10d27s23
              end do                                                     10d27s23
             end do                                                      10d27s23
            else                                                         10d27s23
             do i3=0,irefo(isb)-1                                        10d27s23
              jtmpd=itmpd+i4+irefo(jsb)*i3                              10d27s23
              iadd=id4o(is)+nrowd*(i3+irefo(isb)*i4)                     10d27s23
              do i12=0,nrowd-1                                           10d27s23
               bc(jtmpd)=bc(iadd+i12)                                    10d27s23
               jtmpd=jtmpd+nn                                           10d27s23
              end do                                                     10d27s23
             end do                                                      10d27s23
            end if                                                       10d27s23
           end do                                                        10d27s23
           imap=ibcoff                                                  10d27s23
           ibcoff=imap+nrowd                                            10d27s23
           call enough('dorbd4o.imap',bc,ibc)                           10d27s23
           jmap=imap                                                    10d27s23
           do i4=0,irefo(isb2)-1                                        10d27s23
            i4p=i4+idoubo(isb2)                                         10d27s23
            i3top=i4+isw*(irefo(isb1)-1-i4)                             10d27s23
            do i3=0,i3top                                               10d27s23
             i3p=i3+idoubo(isb1)                                        10d27s23
             irec=i3p+noc(isb1)*i4p                                     10d27s23
             itri=((i4p*(i4p+1))/2)+i3p                                 10d27s23
             itri=itri+isw*(irec-itri)                                  10d27s23
             ibc(jmap)=itri                                             10d27s23
             jmap=jmap+1                                                10d27s23
            end do                                                      10d27s23
           end do                                                       10d27s23
           if(nhere.gt.0)then                                             10d27s23
            itmpi=ibcoff                                                   10d27s23
            ibcoff=itmpi+nmul*noc(jsb)                                   10d27s23
            call enough('dorbd4o.tmpi0',bc,ibc)                          10d27s23
            do iz=itmpi,ibcoff-1                                           10d27s23
             bc(iz)=0d0                                                    10d27s23
            end do                                                         10d27s23
            i0=ioooo(is)                                                   10d27s23
            i10=i1s                                                        10d27s23
            i1n=noc(isb)                                                   10d27s23
            do i2=i2s,i2e                                                  10d27s23
             if(i2.eq.i2e)i1n=i1e                                          10d27s23
             do i1=i10,i1n                                                 10d27s23
              if(i1.gt.idoubo(isb))then                                    10d27s23
               i1m=i1-idoubo(isb)-1                                        10d27s23
               jtmpi=itmpi+i1m+nmul*(i2-1)                                 10d27s23
               do i34=0,nrowd-1                                         10d27s23
                itri=ibc(imap+i34)+i0                                   10d27s23
                bc(jtmpi)=bc(itri)                                        10d27s23
                jtmpi=jtmpi+irefo(isb)                                    10d27s23
               end do                                                   10d27s23
              end if                                                       10d27s23
              i0=i0+nrow                                                   10d27s23
             end do                                                        10d27s23
             i10=1                                                         10d27s23
            end do                                                         10d27s23
            jtt=itt(jsb)+idoubo(jsb)                                       10d27s23
            call dgemm('n','n',irefo(jsb),noc(jsb),nmul,4d0,               10d27s23
     $        bc(itmpd),irefo(jsb),bc(itmpi),nmul,1d0,                  10d27s23
     $        bc(jtt),nbasdws(jsb),'dorbd4o.tt4o')                      10d27s23
            ibcoff=itmpi                                                   10d27s23
           end if                                                        10d27s23
           if(mhere.gt.0)then                                           10d27s23
            mhere2=j2e+1-j2s                                            10d27s23
            itmpi=ibcoff                                                   10d27s23
            ibcoff=itmpi+nmul*mhere2                                    10d27s23
            call enough('dorbd4o.tmpi1',bc,ibc)                          10d27s23
            do iz=itmpi,ibcoff-1                                           10d27s23
             bc(iz)=0d0                                                    10d27s23
            end do                                                         10d27s23
            i0=ionex(is1)                                               10d27s23
            i10=j1s                                                        10d27s23
            i1n=noc(isb)                                                   10d27s23
            do i2=j2s,j2e                                                  10d27s23
             if(i2.eq.j2e)i1n=j1e                                          10d27s23
             do i1=i10,i1n                                                 10d27s23
              if(i1.gt.idoubo(isb))then                                    10d27s23
               i1m=i1-idoubo(isb)-1                                        10d27s23
               jtmpi=itmpi+i1m+nmul*(i2-j2s)                            10d27s23
               do i34=0,nrowd-1                                         10d27s23
                itri=ibc(imap+i34)+i0                                   10d27s23
                bc(jtmpi)=bc(itri)                                        10d27s23
                jtmpi=jtmpi+irefo(isb)                                    10d27s23
               end do                                                   10d27s23
              end if                                                       10d27s23
              i0=i0+nrow                                                   10d27s23
             end do                                                        10d27s23
             i10=1                                                         10d27s23
            end do                                                         10d27s23
            jtt=itt(jsb)+idoubo(jsb)+nbasdws(jsb)*(noc(jsb)+j2s-1)      10d27s23
            call dgemm('n','n',irefo(jsb),mhere2,nmul,4d0,              10d27s23
     $        bc(itmpd),irefo(jsb),bc(itmpi),nmul,1d0,                  10d27s23
     $        bc(jtt),nbasdws(jsb),'dorbd4o.tt1x')                      10d27s23
            ibcoff=itmpi                                                   10d27s23
           end if                                                        10d27s23
           ibcoff=itmpd                                                 10d27s23
          end if                                                        10d27s23
          go to 22
         end if                                                          10d27s23
        end do                                                           10d27s23
   22   continue                                                         10d27s23
       end if                                                           10d27s23
      end do                                                            10d27s23
      return                                                            7d21s23
      end                                                               7d21s23
