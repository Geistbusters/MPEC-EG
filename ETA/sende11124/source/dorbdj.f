c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dorbdj(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,       6d24s24
     $     nsdlk1,isblk1,jmden,ionex,jmats,i3x,bc,ibc,itt,nbasdws)      6d14s24
      implicit real*8 (a-h,o-z)                                         7d31s23
      include "common.store"                                            7d31s23
      dimension idoubo(*),irefo(*),noc(*),nvirt(*),                     6d24s24
     $     isblk(4,*),isblk1(4,*),jmden(*),ionex(*),itt(*),nbasdws(*),  11d1s23
     $     jmats(*),i3x(*),jms(512),ms(28),multh(8,8)                   6d14s24
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data ms/1,5,1,5,1,2,3,6,7,5,5,1,2,3,6,7,1,5,2,3,1,4,5,8,1,6,7,5/  8d1s23
      data multh/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,
     $     4,3,2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,
     $     7,8,5,6,3,4,1,2,8,7,6,5,4,3,2,1/
      ibcoffo=ibcoff                                                    7d31s23
      do is=1,nsdlk                                                     7d31s23
       if(isblk(1,is).eq.isblk(2,is))then                               7d31s23
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              7d31s23
       else                                                             7d31s23
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      7d31s23
       end if                                                           7d31s23
       ncol=nvirt(isblk(3,is))*nvirt(isblk(4,is))                       7d31s23
       jms(is)=ibcoff                                                   7d31s23
       ibcoff=jms(is)+nrow*ncol                                         7d31s23
       call enough('dorbdj.jms',bc,ibc)                                 7d31s23
       do iz=jms(is),ibcoff-1                                           7d31s23
        bc(iz)=0d0                                                      7d31s23
       end do                                                           7d31s23
       call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,      7d31s23
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d31s23
       nhere=ih+1-il                                                    7d31s23
       if(nhere.gt.0)then                                               7d31s23
        i10=i1s                                                         7d31s23
        i1n=nvirt(isblk(3,is))                                          7d31s23
        jj=jmden(is)                                                    7d31s23
        do i2=i2s,i2e                                                   7d31s23
         i2m=i2-1                                                       7d31s23
         if(i2.eq.i2e)i1n=i1e                                           7d31s23
         do i1=i10,i1n                                                  7d31s23
          iadj=jms(is)+nrow*(i1-1+nvirt(isblk(3,is))*i2m)               7d31s23
          jjj=jj                                                        7d31s23
          do i34=0,nrow-1                                               7d31s23
           bc(iadj+i34)=bc(jjj)                                          7d31s23
           jjj=jjj+nhere                                                7d31s23
          end do                                                        7d31s23
          jj=jj+1                                                       7d31s23
         end do                                                         7d31s23
         i10=1                                                          7d31s23
        end do                                                          7d31s23
       end if                                                           7d31s23
      end do                                                            7d31s23
      nwds=ibcoff-ibcoffo                                               7d31s23
      call dws_gsumf(bc(ibcoffo),nwds)                                  7d31s23
      ihit=ibcoff                                                       7d31s23
      ibcoff=ihit+nsdlk                                                 7d31s23
      jhit=ihit-1                                                       7d31s23
      call enough('dorbdj.hit',bc,ibc)                                  7d31s23
      do is=1,nsdlk                                                     7d31s23
       ibc(jhit+is)=0                                                   7d31s23
      end do                                                            7d31s23
      do is=1,nsdlk                                                     7d31s23
       if(ibc(jhit+is).eq.0)then                                        7d31s23
        ncol=nvirt(isblk(3,is))*nvirt(isblk(4,is))                      7d31s23
        if(isblk(1,is).eq.isblk(2,is))then                              7d31s23
         nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2             7d31s23
         do i2=0,nvirt(isblk(4,is))-1                                   7d31s23
          do i1=0,i2-1                                                  7d31s23
           icol12=i1+nvirt(isblk(3,is))*i2                              7d31s23
           iad12=jms(is)+nrow*icol12                                    7d31s23
           icol21=i2+nvirt(isblk(3,is))*i1                              7d31s23
           iad21=jms(is)+nrow*icol21
           do i34=0,nrow-1                                              7d31s23
            avg=0.5d0*(bc(iad12+i34)+bc(iad21+i34))                     7d31s23
            bc(iad12+i34)=avg                                           7d31s23
            bc(iad21+i34)=avg                                           7d31s23
           end do                                                       7d31s23
          end do                                                        7d31s23
         end do                                                         7d31s23
        else                                                            7d31s23
         do js=is,nsdlk                                                 7d31s23
          if(isblk(1,is).eq.isblk(2,js).and.isblk(2,is).eq.isblk(1,js)  7d31s23
     $         .and.isblk(3,is).eq.isblk(4,js))then                     7d31s23
           ibc(jhit+js)=1                                               7d31s23
           do i2=0,nvirt(isblk(4,is))-1                                 7d31s23
            do i1=0,nvirt(isblk(3,is))-1                                7d31s23
             do i4=0,irefo(isblk(2,is))-1                               7d31s23
              do i3=0,irefo(isblk(1,is))-1                              7d31s23
               iad1=jms(is)+i3+irefo(isblk(1,is))*(i4+irefo(isblk(2,is))7d31s23
     $              *(i1+nvirt(isblk(3,is))*i2))                        7d31s23
               iad2=jms(js)+i4+irefo(isblk(2,is))*(i3+irefo(isblk(1,is))7d31s23
     $              *(i2+nvirt(isblk(4,is))*i1))                        7d31s23
               avg=0.5d0*(bc(iad1)+bc(iad2))                            7d31s23
               bc(iad1)=avg                                             7d31s23
               bc(iad2)=avg                                             7d31s23
              end do                                                    7d31s23
             end do                                                     7d31s23
            end do                                                      7d31s23
           end do                                                       7d31s23
           go to 7575                                                   7d31s23
          end if                                                        7d31s23
         end do                                                         7d31s23
        end if                                                          7d31s23
 7575   continue                                                        7d31s23
        ibc(jhit+is)=1                                                  7d31s23
       end if                                                           7d31s23
      end do                                                            7d31s23
c
      do is=1,nsdlk                                                     7d31s23
       if(isblk(1,is).eq.isblk(2,is))then                               7d31s23
        nrow1=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2             7d31s23
        nrow2=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2                 7d31s23
        nrow3=(nvirt(isblk(1,is))*(nvirt(isblk(1,is))+1))/2             11d2s23
        iswitch=0                                                       7d31s23
       else                                                             7d31s23
        nrow1=irefo(isblk(1,is))*irefo(isblk(2,is))                     7d31s23
        nrow2=noc(isblk(1,is))*noc(isblk(2,is))                         7d31s23
        nrow3=nvirt(isblk(1,is))*nvirt(isblk(2,is))                     11d2s23
        iswitch=1                                                       7d31s23
       end if                                                           7d31s23
       ncol=nvirt(isblk(3,is))*nvirt(isblk(4,is))
       imap=ibcoff
       ibcoff=imap+nrow1                                                8d1s23
       call enough('dorbdj.imap',bc,ibc)                                7d31s23
       if(nrow1.gt.0)then                                               7d31s23
        jmap=imap                                                        7d31s23
        do i2=0,irefo(isblk(2,is))-1                                     7d31s23
         i2p=i2+idoubo(isblk(2,is))                                      7d31s23
         itop=i2+iswitch*(irefo(isblk(1,is))-1-i2)                       8d1s23
         do i1=0,itop                                                    7d31s23
          i1p=i1+idoubo(isblk(1,is))                                     7d31s23
          irec=i1p+noc(isblk(1,is))*i2p                                  7d31s23
          itri=((i2p*(i2p+1))/2)+i1p                                     7d31s23
          ibc(jmap+i1)=itri+iswitch*(irec-itri)                          7d31s23
         end do                                                          7d31s23
         jmap=jmap+itop+1                                                7d31s23
        end do                                                           7d31s23
        call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,     7d31s23
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d31s23
        nhere=ih+1-il                                                   7d31s23
        if(nhere.gt.0)then                                              7d31s23
         isb=isblk(1,is)                                                7d31s23
         jsb=isblk(2,is)                                                11d1s23
         if(noc(isb).gt.0)then                                          11d1s23
          nmul=irefo(jsb)*nhere                                         11d1s23
          itmpd=ibcoff                                                  11d1s23
          itmpi=itmpd+irefo(isb)*nmul                                   11d1s23
          ibcoff=itmpi+nmul*noc(isb)                                    11d1s23
          call enough('dorbdj.tmpja',bc,ibc)                            11d1s23
          do iz=itmpi,ibcoff-1                                          11d1s23
           bc(iz)=0d0                                                   11d1s23
          end do                                                        11d1s23
          do i4=0,irefo(jsb)-1                                          11d1s23
           do i3=0,irefo(isb)-1                                         11d1s23
            irec=i3+irefo(isb)*i4                                       11d1s23
            ix=max(i3,i4)                                               11d1s23
            in=min(i3,i4)                                               11d1s23
            itri=((ix*(ix+1))/2)+in                                     11d1s23
            itri=itri+iswitch*(irec-itri)                               11d1s23
            jdj=jms(is)+itri+nrow1*(il-1)                               11d1s23
            jtmpd=itmpd+i3+irefo(isb)*nhere*i4                          11d1s23
            ff=1d0                                                      11d1s23
            if(isb.eq.jsb.and.i3.eq.i4)ff=2d0                           11d1s23
            do i=0,nhere-1                                                11d1s23
             bc(jtmpd)=bc(jdj)*ff                                       11d1s23
             jtmpd=jtmpd+irefo(isb)                                     11d1s23
             jdj=jdj+nrow1                                              11d1s23
            end do
           end do                                                       11d1s23
          end do                                                        11d1s23
          do i12=0,nhere-1                                              11d1s23
           jj=jmats(is)+nrow2*i12                                       11d1s23
           do i4=0,irefo(jsb)-1                                         11d1s23
            i4p=i4+idoubo(jsb)                                          11d1s23
            jtmpi=itmpi+i12+nhere*i4                                    11d1s23
            do i3=0,noc(isb)-1                                          11d1s23
             irec=i3+noc(isb)*i4p                                       11d1s23
             ix=max(i3,i4p)                                             11d1s23
             in=min(i3,i4p)                                             11d1s23
             itri=((ix*(ix+1))/2)+in                                    11d1s23
             itri=itri+iswitch*(irec-itri)+jj                           11d1s23
             bc(jtmpi)=bc(itri)                                         11d1s23
             jtmpi=jtmpi+nmul                                           11d1s23
            end do                                                      11d1s23
           end do                                                       11d1s23
          end do                                                        11d1s23
          jtt=itt(isb)+idoubo(isb)                                      11d1s23
          call dgemm('n','n',irefo(isb),noc(isb),nmul,1d0,              11d1s23
     $         bc(itmpd),irefo(isb),bc(itmpi),nmul,1d0,                 11d1s23
     $         bc(jtt),nbasdws(isb),'dorbdj.tmpja')                     11d1s23
          ibcoff=itmpd                                                  11d1s23
         end if                                                         11d1s23
         if(isb.ne.jsb.and.noc(jsb).gt.0)then                           11d1s23
          nmul=irefo(isb)*nhere                                         11d1s23
          itmpd=ibcoff                                                  11d1s23
          itmpi=itmpd+irefo(jsb)*nmul                                   11d1s23
          ibcoff=itmpi+nmul*noc(jsb)                                    11d1s23
          call enough('dorbdj.tmpjb',bc,ibc)                            11d1s23
          do iz=itmpi,ibcoff-1                                          11d1s23
           bc(iz)=0d0                                                   11d1s23
          end do                                                        11d1s23
          do i4=0,irefo(jsb)-1                                          11d1s23
           do i3=0,irefo(isb)-1                                         11d1s23
            irec=i3+irefo(isb)*i4                                       11d1s23
            jdj=jms(is)+irec+nrow1*(il-1)                               11d1s23
            jtmpd=itmpd+i4+irefo(jsb)*nhere*i3                          11d1s23
            do i=0,nhere-1                                                11d1s23
             bc(jtmpd)=bc(jdj)                                          11d1s23
             jtmpd=jtmpd+irefo(jsb)                                     11d1s23
             jdj=jdj+nrow1                                              11d1s23
            end do
           end do                                                       11d1s23
          end do                                                        11d1s23
          do i12=0,nhere-1                                              11d1s23
           jj=jmats(is)+nrow2*i12                                       11d1s23
           do i4=0,noc(jsb)-1                                           11d1s23
            jtmpi=itmpi+i12+nmul*i4                                     11d2s23
            irec=jj+idoubo(isb)+noc(isb)*i4                                       11d1s23
            do i3=0,irefo(isb)-1                                        11d1s23
             bc(jtmpi)=bc(irec+i3)                                      11d1s23
             jtmpi=jtmpi+nhere                                          11d2s23
            end do                                                      11d1s23
           end do                                                       11d1s23
          end do                                                        11d1s23
          jtt=itt(jsb)+idoubo(jsb)                                      11d1s23
          call dgemm('n','n',irefo(jsb),noc(jsb),nmul,1d0,              11d1s23
     $         bc(itmpd),irefo(jsb),bc(itmpi),nmul,1d0,                 11d1s23
     $         bc(jtt),nbasdws(jsb),'dorbdj.tmpjb')                     11d1s23
          ibcoff=itmpd                                                  11d1s23
         end if                                                         11d1s23
         nv=i2e+1-i2s                                                   11d1s23
         jsb=isblk(4,is)                                                11d1s23
         nmul=nrow1*nvirt(isblk(3,is))                                  11d1s23
         itmpd=ibcoff                                                   11d1s23
         itmpi=itmpd+nvirt(jsb)*nmul                                    11d2s23
         ibcoff=itmpi+nmul*nv                                           11d2s23
         call enough('dorbdj.tmpdjc',bc,ibc)                            11d1s23
         do iz=itmpi,ibcoff-1                                           11d1s23
          bc(iz)=0d0                                                    11d1s23
         end do                                                         11d1s23
         nn=nvirt(jsb)*nrow1                                            11d2s23
         do i2=0,nvirt(jsb)-1                                           11d2s23
          do i1=0,nvirt(isblk(3,is))-1                                  11d2s23
           jtmpd=itmpd+i2+nn*i1                                         11d2s23
           jdj=jms(is)+nrow1*(i1+nvirt(isblk(3,is))*i2)                 11d2s23
           do i34=0,nrow1-1                                             11d1s23
            bc(jtmpd)=bc(jdj+i34)                                       11d2s23
            jtmpd=jtmpd+nvirt(jsb)                                      11d2s23
           end do                                                       11d1s23
          end do                                                        11d1s23
          i10=1                                                         11d1s23
         end do                                                         11d1s23
         i10=i1s                                                        11d1s23
         i1n=nvirt(isblk(3,is))                                         11d1s23
         jj=jmats(is)                                                   11d1s23
         do i2=i2s,i2e                                                  11d1s23
          if(i2.eq.i2e)i1n=i1e                                          11d1s23
          do i1=i10,i1n                                                 11d1s23
           jtmpi=itmpi+nrow1*(i1-1+nvirt(isblk(3,is))*(i2-i2s))         11d1s23
           do i34=0,nrow1-1                                             11d1s23
            kk=ibc(imap+i34)                                            11d1s23
            bc(jtmpi+i34)=bc(kk+jj)                                     11d1s23
           end do                                                       11d1s23
           jj=jj+nrow2                                                  11d1s23
          end do                                                        11d1s23
          i10=1                                                         11d1s23
         end do                                                         11d1s23
         jtt=itt(jsb)+noc(jsb)+nbasdws(jsb)*(noc(jsb)+i2s-1)            11d2s23
         call dgemm('n','n',nvirt(jsb),nv,nmul,2d0,                     11d2s23
     $        bc(itmpd),nvirt(jsb),bc(itmpi),nmul,1d0,                  11d2s23
     $        bc(jtt),nbasdws(jsb),'dorbdj.tmpjc')                      11d1s23
         ibcoff=itmpd                                                   11d1s23
         isb=isblk(2,is)                                                8d1s23
        end if                                                          7d31s23
        do is1=1,nsdlk1
         if(isblk(1,is).eq.isblk1(1,is1).and.                           7d31s23
     $        isblk(2,is).eq.isblk1(2,is1)                              7d31s23
     $       .and.isblk(3,is).eq.isblk1(4,is1))then                      7d31s23
          call ilimts(noc(isblk1(3,is1)),nvirt(isblk1(4,is1)),mynprocg,  7d31s23
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           7d31s23
          isb=isblk(4,is)                                               7d31s23
          nv=i2e+1-i2s                                                   7d31s23
          nhere=ih+1-il
          if(nv.gt.0)then                                               7d31s23
           nmul=nv*nrow1                                                11d3s23
           itmpd=ibcoff                                                 11d3s23
           itmpi=itmpd+nvirt(isb)*nmul                                  11d3s23
           ibcoff=itmpi+nmul*noc(isb)                                   11d3s23
           call enough('dorbdj.tmpd1',bc,ibc)                           11d3s23
           do iz=itmpi,ibcoff-1                                         11d3s23
            bc(iz)=0d0                                                  11d3s23
           end do                                                       11d3s23
           nn=nvirt(isb)*nrow1                                          11d3s23
           do ivp=0,nvirt(isb)-1                                        11d3s23
            do iv=i2s,i2e                                               11d3s23
             jm=jms(is)+nrow1*(iv-1+nvirt(isblk(3,is))*ivp)             11d3s23
             jtmpd=itmpd+ivp+nn*(iv-i2s)                                11d3s23
             do i34=0,nrow1-1                                           11d3s23
              bc(jtmpd)=bc(jm+i34)                                      11d3s23
              jtmpd=jtmpd+nvirt(isb)                                    11d3s23
             end do                                                     11d3s23
            end do                                                      11d3s23
           end do                                                       11d3s23
           i10=i1s                                                      11d3s23
           i1n=noc(isblk1(3,is1))                                       11d3s23
           ii=ionex(is1)                                                11d3s23
           do i2=i2s,i2e                                                11d3s23
            if(i2.eq.i2e)i1n=i1e                                        11d3s23
            do i1=i10,i1n                                               11d3s23
             jtmpi=itmpi+nrow1*(i2-i2s+nv*(i1-1))                       11d3s23
             do i34=0,nrow1-1                                           11d3s23
              kk=ibc(imap+i34)                                          11d3s23
              bc(jtmpi+i34)=bc(ii+kk)                                   11d3s23
             end do                                                     11d3s23
             ii=ii+nrow2                                                11d3s23
            end do                                                      11d3s23
            i10=1                                                       11d3s23
           end do                                                       11d3s23
           jtt=itt(isb)+noc(isb)                                        11d3s23
           call dgemm('n','n',nvirt(isb),noc(isb),nmul,2d0,             11d3s23
     $          bc(itmpd),nvirt(isb),bc(itmpi),nmul,1d0,                11d3s23
     $          bc(jtt),nbasdws(isb),'dorbdj.onex')                     11d3s23
           ibcoff=itmpd                                                 11d3s23
          end if                                                        7d31s23
         end if                                                          7d31s23
        end do                                                           7d31s23
       end if                                                           7d31s23
      end do                                                            7d31s23
c
c     for 3virt integrals, symmetrize virt indicies in densities.
c
      ngt=0                                                             8d1s23
      nlt=nsdlk+1                                                       8d1s23
      do is=1,nsdlk                                                     8d1s23
       if(isblk(3,is).ge.isblk(4,is))then                               8d1s23
        ngt=max(ngt,is)                                                 8d1s23
       else                                                             8d1s23
        nlt=min(nlt,is)                                                 8d1s23
       end if                                                           8d1s23
      end do                                                            8d1s23
      if(ngt.gt.nlt)then                                                8d1s23
       write(6,*)('ooops!')                                              8d1s23
       stop 'dorbdj'                                                    8d1s23
      end if                                                            8d1s23
      do is=1,ngt                                                       8d1s23
       if(isblk(3,is).ne.isblk(4,is))then                               8d1s23
        nrow=irefo(isblk(1,is))*irefo(isblk(2,is))                      8d1s23
        do js=nlt,nsdlk                                                 8d1s23
         if(isblk(1,is).eq.isblk(1,js).and.isblk(2,is).eq.isblk(2,js)   8d1s23
     $        .and.isblk(3,is).eq.isblk(4,js))then                      8d1s23
          ipart=js                                                      8d1s23
          go to 888                                                     8d1s23
         end if                                                         8d1s23
        end do                                                          8d1s23
        write(6,*)('could not find partner for '),(isblk(j,is),j=1,4)   8d1s23
        stop 'dorbdj'                                                   8d1s23
  888   continue                                                        8d1s23
        ncol=nvirt(isblk(3,is))*nvirt(isblk(4,is))                      8d1s23
        do i2=0,nvirt(isblk(4,is))-1                                    8d1s23
         do i1=0,nvirt(isblk(3,is))-1                                   8d1s23
          i12=jms(is)+nrow*(i1+nvirt(isblk(3,is))*i2)                   8d2s23
          i21=jms(ipart)+nrow*(i2+nvirt(isblk(4,is))*i1)                8d2s23
          do i34=0,nrow-1                                               8d1s23
           sum=bc(i12+i34)+bc(i21+i34)                                  8d2s23
           bc(i12+i34)=sum                                              8d2s23
          end do                                                        8d1s23
         end do                                                         8d1s23
        end do                                                          8d1s23
        iswitch=1                                                       8d1s23
       else                                                             8d1s23
        nrow=(irefo(isblk(1,is))*(irefo(isblk(1,is))+1))/2              8d1s23
        ncol=(nvirt(isblk(3,is))*(nvirt(isblk(3,is))+1))/2              8d1s23
        itmp=ibcoff                                                     8d1s23
        ibcoff=itmp+nrow*ncol                                           8d1s23
        call enough('dorbdj.tmp.den',bc,ibc)                            8d1s23
        jtmp=itmp                                                       8d1s23
        do i2=0,nvirt(isblk(4,is))-1                                    8d1s23
         i2p=i2+idoubo(isb)+irefo(isb)+1
         do i1=0,i2-1                                                   8d1s23
          i12=jms(is)+nrow*(i1+nvirt(isblk(4,is))*i2)                   8d1s23
          i21=jms(is)+nrow*(i2+nvirt(isblk(4,is))*i1)                   8d1s23
          i1p=i1+idoubo(isb)+irefo(isb)+1
          i34=0
          do i4=0,irefo(isb)-1
           i4p=i4+idoubo(isb)+1
           do i3=0,i4
            i3p=i3+idoubo(isb)+1
            if(max(abs(bc(i12+i34)),abs(bc(i21+i34))).gt.1d-10)then
             icd=multh(ms(i3p),ms(i4p))
             sum=bc(i12+i34)+bc(i21+i34)
            end if
            i34=i34+1
           end do
          end do
          do i34=0,nrow-1                                               8d1s23
           bc(jtmp+i34)=bc(i12+i34)+bc(i21+i34)                         8d1s23
          end do                                                        8d1s23
          jtmp=jtmp+nrow                                                8d1s23
         end do                                                         8d1s23
         i22=jms(is)+nrow*(i2+nvirt(isblk(4,is))*i2)                    8d1s23
         do i34=0,nrow-1                                                8d1s23
          bc(jtmp+i34)=bc(i22+i34)                                      8d1s23
         end do                                                         8d1s23
         jtmp=jtmp+nrow                                                 8d1s23
        end do                                                          8d1s23
        do i=0,nrow*ncol-1                                              8d1s23
         bc(jms(is)+i)=bc(itmp+i)                                       8d1s23
        end do                                                          8d1s23
        ibcoff=itmp                                                     8d1s23
        iswitch=0                                                       8d1s23
       end if                                                           8d1s23
       if(min(nrow,ncol).gt.0)then                                      8d1s23
        do is1=1,nsdlk1                                                  8d1s23
         if(isblk(2,is).eq.isblk1(3,is1)                                8d1s23
     $        .and.isblk(3,is).eq.isblk1(1,is1)                         8d1s23
     $       .and.isblk(4,is).eq.isblk1(2,is1))then                     8d1s23
          call ilimts(noc(isblk1(3,is1)),nvirt(isblk1(4,is1)),mynprocg,  8d1s23
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d1s23
          nhere=ih+1-il                                                  8d1s23
          if(nhere.gt.0)then                                             8d1s23
           isb=isblk(1,is)                                               8d1s23
           nv=i2e+1-i2s                                                  8d1s23
           if(irefo(isblk1(3,is1)).gt.0)then                            11d2s23
            nmul=ncol*irefo(isblk1(3,is1))                              11d2s23
            itmpd=ibcoff                                                 11d2s23
            itmpi=itmpd+irefo(isb)*nmul                                 11d2s23
            ibcoff=itmpi+nmul*nv                                        11d2s23
            call enough('dorbdj.3xa',bc,ibc)                            11d2s23
            do iz=itmpi,ibcoff-1                                        11d2s23
             bc(iz)=0d0                                                 11d2s23
            end do                                                      11d2s23
            nn=irefo(isb)*ncol                                          11d2s23
            do i4=0,irefo(isblk(2,is))-1                                11d2s23
             do i3=0,irefo(isb)-1                                       11d2s23
              irec=i3+irefo(isblk(1,is))*i4                             11d2s23
              ix=max(i3,i4)                                             11d2s23
              in=min(i3,i4)                                             11d2s23
              itri=((ix*(ix+1))/2)+in                                   11d2s23
              itri=itri+iswitch*(irec-itri)                             11d2s23
              ff=1d0                                                    11d2s23
              if(iswitch.eq.0.and.i3.eq.i4)ff=2d0                       11d3s23
              jtmpd=itmpd+i3+nn*i4                                      11d2s23
              ifrm=jms(is)+itri                                         11d2s23
              do i12=0,ncol-1                                           11d2s23
               bc(jtmpd)=bc(ifrm)*ff                                    11d2s23
               jtmpd=jtmpd+irefo(isb)                                   11d2s23
               ifrm=ifrm+nrow                                           11d2s23
              end do                                                    11d2s23
             end do                                                     11d2s23
            end do                                                      11d2s23
            i10=i1s                                                     11d2s23
            i1n=noc(isblk1(3,is1))                                      11d2s23
            j3x=i3x(is1)                                                11d2s23
            do i2=i2s,i2e                                               11d2s23
             if(i2.eq.i2e)i1n=i1e                                       11d2s23
             do i1=i10,i1n                                              11d2s23
              if(i1.gt.idoubo(isblk1(3,is1)))then                       11d2s23
               jtmpi=itmpi+ncol*(i1-1-idoubo(isblk1(3,is1))
     $             +irefo(isblk1(3,is1))*(i2-i2s))                      11d2s23
               do i34=0,ncol-1                                           11d2s23
                bc(jtmpi+i34)=bc(j3x+i34)                               11d2s23
               end do                                                    11d2s23
              end if                                                    11d2s23
              j3x=j3x+ncol                                              11d2s23
             end do                                                     11d2s23
             i10=1                                                      11d2s23
            end do                                                      11d2s23
            jtt=itt(isb)+idoubo(isb)+nbasdws(isb)*(noc(isb)+i2s-1)      11d2s23
            call dgemm('n','n',irefo(isb),nv,nmul,1d0,                  11d2s23
     $           bc(itmpd),irefo(isb),bc(itmpi),nmul,1d0,               11d2s23
     $           bc(jtt),nbasdws(isb),'dorbdj.3xa')                     11d2s23
            ibcoff=itmpd                                                11d2s23
           end if                                                       11d2s23
          end if                                                         8d1s23
         else if(isblk(2,is).eq.isblk1(4,is1)                                8d1s23
     $        .and.isblk(3,is).eq.isblk1(1,is1)                         8d1s23
     $       .and.isblk(4,is).eq.isblk1(2,is1).and.iswitch.ne.0)then    11d3s23
          call ilimts(noc(isblk1(3,is1)),nvirt(isblk1(4,is1)),mynprocg,  8d1s23
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          8d1s23
          nhere=ih+1-il                                                  8d1s23
          if(nhere.gt.0)then                                             8d1s23
           isb=isblk(2,is)                                               8d1s23
           nv=i2e+1-i2s                                                  8d1s23
           if(irefo(isblk1(3,is1)).gt.0)then                            11d2s23
            nmul=ncol*irefo(isblk1(3,is1))                              11d2s23
            itmpd=ibcoff                                                 11d2s23
            itmpi=itmpd+irefo(isb)*nmul                                 11d2s23
            ibcoff=itmpi+nmul*nv                                        11d2s23
            call enough('dorbdj.3xb',bc,ibc)                            11d2s23
            do iz=itmpi,ibcoff-1                                        11d2s23
             bc(iz)=0d0                                                 11d2s23
            end do                                                      11d2s23
            nn=irefo(isb)*ncol                                          11d2s23
            do i4=0,irefo(isb)-1                                        11d3s23
             do i3=0,irefo(isblk(1,is))-1                               11d3s23
              irec=i3+irefo(isblk(1,is))*i4                             11d2s23
              jtmpd=itmpd+i4+nn*i3                                      11d2s23
              ifrm=jms(is)+irec                                         11d3s23
              do i12=0,ncol-1                                           11d2s23
               bc(jtmpd)=bc(ifrm)                                       11d3s23
               jtmpd=jtmpd+irefo(isb)                                   11d2s23
               ifrm=ifrm+nrow                                           11d2s23
              end do                                                    11d2s23
             end do                                                     11d2s23
            end do                                                      11d2s23
            i10=i1s                                                     11d2s23
            i1n=noc(isblk1(3,is1))                                      11d2s23
            j3x=i3x(is1)                                                11d2s23
            do i2=i2s,i2e                                               11d2s23
             if(i2.eq.i2e)i1n=i1e                                       11d2s23
             do i1=i10,i1n                                              11d2s23
              if(i1.gt.idoubo(isblk1(3,is1)))then                       11d2s23
               jtmpi=itmpi+ncol*(i1-1-idoubo(isblk1(3,is1))
     $             +irefo(isblk1(3,is1))*(i2-i2s))                      11d2s23
               do i34=0,ncol-1                                           11d2s23
                bc(jtmpi+i34)=bc(j3x+i34)                               11d2s23
               end do                                                    11d2s23
              end if                                                    11d2s23
              j3x=j3x+ncol                                              11d2s23
             end do                                                     11d2s23
             i10=1                                                      11d2s23
            end do                                                      11d2s23
            jtt=itt(isb)+idoubo(isb)+nbasdws(isb)*(noc(isb)+i2s-1)      11d2s23
            call dgemm('n','n',irefo(isb),nv,nmul,1d0,                  11d2s23
     $           bc(itmpd),irefo(isb),bc(itmpi),nmul,1d0,               11d2s23
     $           bc(jtt),nbasdws(isb),'dorbdj.3xb')                     11d2s23
            ibcoff=itmpd                                                11d2s23
           end if                                                       11d2s23
          end if
         end if                                                          8d1s23
        end do                                                           8d1s23
       end if                                                           8d1s23
      end do                                                            8d1s23
      ibcoff=ibcoffo                                                    7d31s23
      return
      end
