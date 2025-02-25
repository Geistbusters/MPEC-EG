c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dorbdk(nsymb,idoubo,irefo,noc,nvirt,nsdlkk,            6d24s24
     $     isblkk,nsdlk1,isblk1,kmden,ionex,kmats,i3x,bc,ibc,itt,       6d14s24
     $     nbasdws)                                                     11d6s23
      implicit real*8 (a-h,o-z)                                         7d31s23
      include "common.store"                                            7d31s23
      dimension idoubo(*),irefo(*),noc(*),nvirt(*),                     6d24s24
     $     isblkk(4,*),isblk1(4,*),kmden(*),ionex(*),itt(*),nbasdws(*), 11d6s23
     $     kmats(*),i3x(*),kms(512),ms(28),multh(8,8)                   6d14s24
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      data ms/1,5,1,5,1,2,3,6,7,5,5,1,2,3,6,7,1,5,2,3,1,4,5,8,1,6,7,5/  8d1s23
      data multh/1,2,3,4,5,6,7,8,2,1,4,3,6,5,8,7,3,4,1,2,7,8,5,6,
     $     4,3,2,1,8,7,6,5,5,6,7,8,1,2,3,4,6,5,8,7,2,1,4,3,
     $     7,8,5,6,3,4,1,2,8,7,6,5,4,3,2,1/
      ibcoffo=ibcoff                                                    7d31s23
      do is=1,nsdlkk                                                     7d31s23
       nrow=irefo(isblkk(1,is))*irefo(isblkk(2,is))                      7d31s23
       ncol=nvirt(isblkk(3,is))*nvirt(isblkk(4,is))                       7d31s23
       kms(is)=ibcoff                                                   7d31s23
       ibcoff=kms(is)+nrow*ncol                                         7d31s23
       call enough('dorbdk.kms',bc,ibc)                                 7d31s23
       do iz=kms(is),ibcoff-1                                           7d31s23
        bc(iz)=0d0                                                      7d31s23
       end do                                                           7d31s23
       call ilimts(nvirt(isblkk(3,is)),nvirt(isblkk(4,is)),mynprocg,      7d31s23
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            7d31s23
       nhere=ih+1-il                                                    7d31s23
       if(nhere.gt.0)then                                               7d31s23
        i10=i1s                                                         7d31s23
        i1n=nvirt(isblkk(3,is))                                          7d31s23
        jj=kmden(is)                                                    7d31s23
        do i2=i2s,i2e                                                   7d31s23
         i2m=i2-1                                                       7d31s23
         if(i2.eq.i2e)i1n=i1e                                           7d31s23
         do i1=i10,i1n                                                  7d31s23
          iadj=kms(is)+nrow*(i1-1+nvirt(isblkk(3,is))*i2m)               7d31s23
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
      ihit=ibcoff                                                       8d2s23
      ibcoff=ihit+nsdlkk                                                8d2s23
      call enough('dorbdk.hit',bc,ibc)                                  8d2s23
      jhit=ihit-1                                                       8d2s23
      do i=1,nsdlkk                                                     8d2s23
       ibc(jhit+i)=0                                                    8d2s23
      end do                                                            8d2s23
      do is=1,nsdlkk                                                    8d2s23
       nrow=irefo(isblkk(1,is))*irefo(isblkk(2,is))                     8d2s23
       nrow2=noc(isblkk(1,is))*noc(isblkk(2,is))                        8d2s23
       ncol=nvirt(isblkk(3,is))*nvirt(isblkk(4,is))                     8d2s23
       if(ibc(jhit+is).eq.0)then                                        8d2s23
        if(isblkk(1,is).eq.isblkk(2,is))then                            8d2s23
         do i2=0,nvirt(isblkk(4,is))-1                                  8d2s23
          do i1=0,nvirt(isblkk(3,is))-1                                 8d2s23
           do i4=0,irefo(isblkk(2,is))-1                                8d2s23
            do i3=0,irefo(isblkk(1,is))-1                               8d2s23
             ij=kms(is)+i3+irefo(isblkk(1,is))*(i4+irefo(isblkk(2,is))  8d2s23
     $            *(i1+nvirt(isblkk(3,is))*i2))                         8d2s23
             ji=kms(is)+i4+irefo(isblkk(1,is))*(i3+irefo(isblkk(2,is))  8d2s23
     $            *(i2+nvirt(isblkk(3,is))*i1))                         8d2s23
             avg=0.5d0*(bc(ij)+bc(ji))                                  8d2s23
             bc(ij)=avg                                                 8d2s23
             bc(ji)=avg                                                 8d2s23
            end do                                                      8d2s23
           end do                                                       8d2s23
          end do                                                        8d2s23
         end do                                                         8d2s23
        else                                                            8d2s23
         do js=is+1,nsdlkk                                              8d2s23
          if(isblkk(1,is).eq.isblkk(2,js).and.                          8d2s23
     $       isblkk(2,is).eq.isblkk(1,js).and.                          8d2s23
     $       isblkk(3,is).eq.isblkk(4,js))then                          8d2s23
           ibc(jhit+js)=1                                               8d2s23
           do i2=0,nvirt(isblkk(4,is))-1                                  8d2s23
            do i1=0,nvirt(isblkk(3,is))-1                                 8d2s23
             do i4=0,irefo(isblkk(2,is))-1                                8d2s23
              do i3=0,irefo(isblkk(1,is))-1                               8d2s23
               ij=kms(is)+i3+irefo(isblkk(1,is))*(i4+irefo(isblkk(2,is))  8d2s23
     $            *(i1+nvirt(isblkk(3,is))*i2))                         8d2s23
               ji=kms(js)+i4+irefo(isblkk(2,is))*(i3+irefo(isblkk(1,is))  8d2s23
     $            *(i2+nvirt(isblkk(4,is))*i1))                         8d2s23
               avg=0.5d0*(bc(ij)+bc(ji))                                  8d2s23
   22          format(4i4,3es15.7)
               bc(ij)=avg                                                 8d2s23
               bc(ji)=avg                                                 8d2s23
              end do                                                      8d2s23
             end do                                                       8d2s23
            end do                                                        8d2s23
           end do                                                         8d2s23
          end if                                                        8d2s23
         end do                                                         8d2s23
        end if                                                          8d2s23
       end if                                                           8d2s23
       call ilimts(nvirt(isblkk(3,is)),nvirt(isblkk(4,is)),mynprocg,     8d2s23
     $      mynowprog,il,ih,i1s,i1e,i2s,i2e)                            8d2s23
       nhere=ih+1-il                                                    8d2s23
       isb=isblkk(1,is)                                                 11d7s23
       if(min(nhere,irefo(isb),irefo(isblkk(2,is))).gt.0)then           11d7s23
        nmul=nhere*irefo(isblkk(2,is))                                  11d7s23
        itmpd=ibcoff                                                    11d7s23
        itmpi=itmpd+irefo(isb)*nmul                                     11d7s23
        ibcoff=itmpi+nmul*noc(isb)                                      11d7s23
        call enough('dorbdk.ka',bc,ibc)                                 11d7s23
        do iz=itmpi,ibcoff-1                                            11d7s23
         bc(iz)=0d0                                                     11d7s23
        end do                                                          11d7s23
        do i12=0,nhere-1                                                11d7s23
         do i4=0,irefo(isblkk(2,is))-1                                  11d7s23
          jk=kms(is)+irefo(isb)*(i4+irefo(isblkk(2,is))*(il+i12-1))     11d7s23
          jtmpd=itmpd+irefo(isb)*(i4+irefo(isblkk(2,is))*i12)           11d7s23
          do i3=0,irefo(isb)-1                                          11d7s23
           bc(jtmpd+i3)=bc(jk+i3)                                       11d7s23
          end do                                                        11d7s23
         end do                                                         11d7s23
        end do                                                          11d7s23
        i10=i1s                                                         11d7s23
        i1n=nvirt(isblkk(3,is))                                         11d7s23
        kk=kmats(is)                                                    11d7s23
        icol=0                                                          11d7s23
        do i2=i2s,i2e                                                   11d7s23
         if(i2.eq.i2e)i1n=i1e                                           11d7s23
         do i1=i10,i1n                                                  11d7s23
          do i4=0,irefo(isblkk(2,is))-1                                 11d7s23
           i4p=i4+idoubo(isblkk(2,is))                                  11d7s23
           kkk=kk+noc(isb)*i4p                                          11d7s23
           jtmpi=itmpi+i4+irefo(isblkk(2,is))*icol                      11d7s23
           do i3=0,noc(isb)-1                                           11d7s23
            bc(jtmpi)=bc(kkk+i3)                                        11d7s23
            jtmpi=jtmpi+nmul                                            11d7s23
           end do                                                       11d7s23
          end do                                                        11d7s23
          kk=kk+nrow2                                                   11d7s23
          icol=icol+1                                                   11d7s23
         end do                                                         11d7s23
         i10=1                                                          11d7s23
        end do                                                          11d7s23
        jtt=itt(isb)+idoubo(isb)                                        11d7s23
        call dgemm('n','n',irefo(isb),noc(isb),nmul,2d0,                11d7s23
     $       bc(itmpd),irefo(isb),bc(itmpi),nmul,1d0,                   11d7s23
     $       bc(jtt),nbasdws(isb),'dorbdk.ka')                          11d7s23
        ibcoff=itmpd                                                    6d20s24
        nv=i2e+1-i2s                                                    11d7s23
        if(nv.gt.0)then
         nmul=nrow*nvirt(isblkk(3,is))                                  11d7s23
         itmpd=ibcoff                                                   11d7s23
         itmpi=itmpd+nvirt(isblkk(4,is))*nmul                           11d7s23
         ibcoff=itmpi+nmul*nv                                           11d7s23
         call enough('dorbdk.kb',bc,ibc)                                11d7s23
         do iz=itmpi,ibcoff-1                                           11d7s23
          bc(iz)=0d0                                                    11d7s23
         end do                                                         11d7s23
         nn=nvirt(isblkk(4,is))*nrow                                    11d7s23
         do ivp=0,nvirt(isblkk(4,is))-1                                 11d7s23
          do iv=0,nvirt(isblkk(3,is))-1                                 11d7s23
           kd=kms(is)+nrow*(iv+nvirt(isblkk(3,is))*ivp)                 11d7s23
           jtmpd=itmpd+ivp+nn*iv                                        11d7s23
           do i12=0,nrow-1                                              11d7s23
            bc(jtmpd)=bc(kd+i12)                                        11d7s23
            jtmpd=jtmpd+nvirt(isblkk(4,is))                             11d7s23
           end do                                                       11d7s23
          end do                                                        11d7s23
         end do                                                         11d7s23
         i10=i1s                                                        11d7s23
         i1n=nvirt(isblkk(3,is))                                        11d7s23
         kk=kmats(is)                                                   11d7s23
         do i2=i2s,i2e                                                  11d7s23
          if(i2.eq.i2e)i1n=i1e                                          11d7s23
          do i1=i10,i1n                                                 11d7s23
           jtmpi=itmpi+nrow*(i1-1+nvirt(isblkk(3,is))*(i2-i2s))         11d7s23
           do i4=0,irefo(isblkk(2,is))-1                                11d7s23
            i4p=i4+idoubo(isblkk(2,is))                                 11d7s23
            kkk=kk+idoubo(isblkk(1,is))+noc(isblkk(1,is))*i4p           11d7s23
            do i3=0,irefo(isblkk(1,is))-1                               11d7s23
             bc(jtmpi+i3)=bc(kkk+i3)                                    11d7s23
            end do                                                      11d7s23
            jtmpi=jtmpi+irefo(isblkk(1,is))                             11d7s23
           end do                                                       11d7s23
           kk=kk+nrow2                                                  11d7s23
          end do                                                        11d7s23
          i10=1                                                         11d7s23
         end do                                                         11d7s23
         jtt=itt(isblkk(4,is))+noc(isblkk(4,is))+nbasdws(isblkk(4,is))* 11d7s23
     $        (noc(isblkk(4,is))+i2s-1)                                 11d7s23
         call dgemm('n','n',nvirt(isblkk(4,is)),nv,nmul,2d0,            11d7s23
     $        bc(itmpd),nvirt(isblkk(4,is)),bc(itmpi),nmul,1d0,         11d7s23
     $        bc(jtt),nbasdws(isblkk(4,is)),'dorbdk.kb')                11d7s23
         ibcoff=itmpd
        end if
       end if                                                           8d2s23
       do is1=1,nsdlk1                                                  8d2s23
        if(isblk1(1,is1).eq.isblk1(2,is1))then                          11d6s23
         nrow1=(noc(isblk1(1,is1))*(noc(isblk1(1,is1))+1))/2            11d6s23
         nrow3=(nvirt(isblk1(1,is1))*(nvirt(isblk1(1,is1))+1))/2        11d6s23
         isw=0                                                          11d6s23
        else                                                            11d6s23
         nrow1=noc(isblk1(1,is1))*noc(isblk1(2,is1))                    11d6s23
         nrow3=nvirt(isblk1(1,is1))*nvirt(isblk1(2,is1))                11d6s23
         isw=1                                                          11d6s23
        end if                                                          11d6s23
        call ilimts(noc(isblk1(3,is1)),nvirt(isblk1(4,is1)),mynprocg,   11d6s23
     $       mynowprog,il,ih,i1s,i1e,i2s,i2e)                           11d6s23
        nhere=ih+1-il                                                   11d6s23
        nv=i2e+1-i2s                                                    11d6s23
        if(min(nhere,nv).gt.0)then                                      11d6s23
         if(isblkk(2,is).eq.isblk1(3,is1).and.                          11d7s23
     $     isblkk(3,is).eq.isblk1(4,is1))then                           11d7s23
          nmul=nrow*nv                                                  11d6s23
          imap=ibcoff                                                    11d6s23
          itmpd=imap+noc(isblkk(4,is))*irefo(isblkk(1,is))              11d7s23
          itmpi=itmpd+nvirt(isblkk(4,is))*nmul                          11d6s23
          ibcoff=itmpi+nmul*noc(isblkk(4,is))                           11d7s23
          call enough('dorbdk.1xa',bc,ibc)                              11d6s23
          jmap=0
          igot=0                                                        11d7s23
          if(isblkk(1,is).eq.isblk1(1,is1))then                          11d6s23
           igot=1                                                       11d7s23
           jmap=imap                                                    11d6s23
           do i4=0,noc(isblk1(2,is1))-1                                  11d6s23
            do i3=0,irefo(isblk1(1,is1))-1                               11d6s23
             i3p=i3+idoubo(isblk1(1,is1))                               11d6s23
             irec=i3p+noc(isblk1(1,is1))*i4                             11d6s23
             ix=max(i4,i3p)                                             11d6s23
             in=min(i4,i3p)                                             11d6s23
             itri=((ix*(ix+1))/2)+in                                    11d6s23
             itri=itri+isw*(irec-itri)                                  11d6s23
             ibc(jmap+i3)=itri                                          11d6s23
            end do                                                      11d6s23
            jmap=jmap+irefo(isblk1(1,is1))                              11d6s23
           end do                                                       11d6s23
          else if(isblkk(1,is).eq.isblk1(2,is1))then                     11d6s23
           igot=1                                                       11d7s23
           jmap=imap                                                    11d6s23
           do i4=0,noc(isblk1(1,is1))-1                                 11d6s23
            do i3=0,irefo(isblk1(2,is1))-1                              11d6s23
             i3p=i3+idoubo(isblk1(2,is1))                               11d6s23
             irec=i4+noc(isblk1(1,is1))*i3p                             11d6s23
             ibc(jmap+i3)=irec                                          11d6s23
            end do                                                      11d6s23
            jmap=jmap+irefo(isblk1(2,is1))                              11d6s23
           end do                                                       11d6s23
          end if                                                         11d6s23
          if(jmap.gt.imap)then                                          11d6s23
           do iz=itmpi,ibcoff-1                                          11d6s23
            bc(iz)=0d0                                                   11d6s23
           end do                                                        11d6s23
           nn=nvirt(isblkk(4,is))*nrow                                  11d7s23
           do i2=i2s,i2e                                                 11d6s23
            do iv=0,nvirt(isblkk(4,is))-1                               11d7s23
             jdk=kms(is)+nrow*(i2-1+nvirt(isblkk(3,is))*iv)             11d7s23
             jtmpd=itmpd+iv+nn*(i2-i2s)                                  11d6s23
             do i34=0,nrow-1                                             11d6s23
              bc(jtmpd)=bc(jdk+i34)                                      11d6s23
              jtmpd=jtmpd+nvirt(isblkk(4,is))                           11d7s23
             end do                                                      11d6s23
            end do                                                       11d6s23
           end do                                                        11d6s23
           i10=i1s                                                       11d6s23
           i1n=noc(isblk1(3,is1))                                        11d6s23
           i1x=ionex(is1)                                                11d6s23
           do i2=i2s,i2e                                                 11d6s23
            if(i2.eq.i2e)i1n=i1e                                         11d6s23
            do i1=i10,i1n                                                11d6s23
             if(i1.gt.idoubo(isblk1(3,is1)))then                        11d7s23
              i1m=i1-1-idoubo(isblk1(3,is1))                             11d6s23
              jmap=imap                                                   11d6s23
              do i4=0,noc(isblkk(4,is))-1                               11d7s23
               jtmpi=itmpi+irefo(isblkk(1,is))*(i1m                      11d6s23
     $              +irefo(isblkk(2,is))*(i2-i2s+nv*i4))                 11d6s23
               do i3=0,irefo(isblkk(1,is))-1                              11d6s23
                kk=ibc(jmap+i3)                                          11d6s23
                bc(jtmpi+i3)=bc(i1x+kk)                                  11d6s23
               end do                                                    11d6s23
               jmap=jmap+irefo(isblkk(1,is))                             11d6s23
              end do                                                     11d6s23
             end if                                                      11d6s23
             i1x=i1x+nrow1                                               11d6s23
            end do                                                       11d6s23
            i10=1                                                        11d6s23
           end do                                                        11d6s23
           jtt=itt(isblkk(4,is))+noc(isblkk(4,is))                       11d6s23
           call dgemm('n','n',nvirt(isblkk(4,is)),noc(isblkk(4,is)),    11d7s23
     $          nmul,2d0,bc(itmpd),nvirt(isblkk(4,is)),bc(itmpi),nmul,  11d7s23
     $          1d0,bc(jtt),nbasdws(isblkk(4,is)),'dorbdk.1xa')         11d7s23
           ibcoff=imap                                                  11d7s23
          end if                                                        11d7s23
          if(min(igot,irefo(isblkk(2,is)),irefo(isblkk(1,is))).gt.0)then11d7s23
           nmul=nvirt(isblkk(4,is))*irefo(isblkk(2,is))*nv              11d7s23
           itmpd=ibcoff                                                 11d7s23
           itmpi=itmpd+irefo(isblkk(1,is))*nmul                         11d7s23
           ibcoff=itmpi+nmul*nvirt(isblkk(1,is))                        11d7s23
           call enough('dorbdk.3x',bc,ibc)                              11d7s23
           do iz=itmpi,ibcoff-1                                         11d7s23
            bc(iz)=0d0                                                  11d7s23
           end do                                                       11d7s23
           do ivp=0,nvirt(isblkk(4,is))-1                               11d7s23
            do i2=i2s,i2e                                               11d7s23
             kk=kms(is)+nrow*(i2-1+nvirt(isblkk(3,is))*ivp)             11d7s23
             jtmpd=itmpd+nrow*(i2-i2s+nv*ivp)                           11d7s23
             do i34=0,nrow-1                                            11d7s23
              bc(jtmpd+i34)=bc(kk+i34)                                  11d7s23
             end do                                                     11d7s23
            end do                                                      11d7s23
           end do                                                       11d7s23
           i10=i1s                                                      11d7s23
           i1n=noc(isblkk(2,is))                                        11d7s23
           j3x=i3x(is1)                                                 11d7s23
           if(isblkk(1,is).eq.isblk1(1,is1))then                        11d7s23
            do i2=i2s,i2e                                                11d7s23
             if(i2.eq.i2e)i1n=i1e                                        11d7s23
             do i1=i10,i1n                                               11d7s23
              if(i1.gt.idoubo(isblk1(3,is1)))then                       11d7s23
               i1m=i1-idoubo(isblk1(3,is1))-1                           11d7s23
               do i4=0,nvirt(isblk1(2,is1))-1                           11d7s23
                jtmpi=itmpi+i1m+irefo(isblk1(3,is1))*(i2-i2s+nv*i4)     11d7s23
                do i3=0,nvirt(isblk1(1,is1))-1                          11d7s23
                 irec=i3+nvirt(isblk1(1,is1))*i4                        11d7s23
                 ix=max(i3,i4)                                          11d7s23
                 in=min(i3,i4)                                          11d7s23
                 itri=((ix*(ix+1))/2)+in                                11d7s23
                 itri=itri+isw*(irec-itri)+j3x                          11d7s23
                 bc(jtmpi)=bc(itri)                                     11d7s23
                 jtmpi=jtmpi+nmul                                       11d7s23
                end do                                                  11d7s23
               end do                                                   11d7s23
              end if                                                    11d7s23
              j3x=j3x+nrow3                                              11d7s23
             end do                                                      11d7s23
             i10=1                                                       11d7s23
            end do
           else if(isblkk(1,is).eq.isblk1(2,is1))then                   11d7s23
            do i2=i2s,i2e                                                11d7s23
             if(i2.eq.i2e)i1n=i1e                                        11d7s23
             do i1=i10,i1n                                               11d7s23
              if(i1.gt.idoubo(isblk1(3,is1)))then                       11d7s23
               i1m=i1-idoubo(isblk1(3,is1))-1                           11d7s23
               do i4=0,nvirt(isblk1(2,is1))-1                           11d7s23
                do i3=0,nvirt(isblk1(1,is1))-1                          11d7s23
                 irec=i3+nvirt(isblk1(1,is1))*i4+j3x                    11d7s23
                 jtmpi=itmpi+i1m+irefo(isblk1(3,is1))*(i2-i2s+nv*(i3    11d7s23
     $                +nvirt(isblk1(1,is1))*i4))                        11d7s23
                 bc(jtmpi)=bc(irec)                                     11d7s23
                end do                                                  11d7s23
               end do                                                   11d7s23
              end if                                                    11d7s23
              j3x=j3x+nrow3                                              11d7s23
             end do                                                      11d7s23
             i10=1                                                       11d7s23
            end do
           end if                                                       11d7s23
           jtt=itt(isblkk(1,is))+idoubo(isblkk(1,is))                   11d7s23
     $          +nbasdws(isblkk(1,is))*noc(isblkk(1,is))                11d7s23
           call dgemm('n','n',irefo(isblkk(1,is)),nvirt(isblkk(1,is)),  11d7s23
     $        nmul,2d0,bc(itmpd),irefo(isblkk(1,is)),bc(itmpi),nmul,1d0,11d7s23
     $          bc(jtt),nbasdws(isblkk(1,is)),'dorbdk.3x')              11d7s23
          end if                                                        11d6s23
          ibcoff=imap                                                   11d6s23
         end if                                                          11d6s23
        end if                                                          11d6s23
       end do                                                           8d2s23
      end do                                                            8d2s23
      ibcoff=ibcoffo
      return
      end
