c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dorbd3x(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,      11d14s23
     $     nsdlk1,isblk1,nsdlkk,isblkk,id3x,jmats,kmats,i3x,itt,nbasdws,11d14s23
     $     bc,ibc)                                                      11d14s23
      implicit real*8 (a-h,o-z)
      dimension idoubo(*),irefo(*),noc(*),nvirt(*),                     11d14s23
     $     isblk(4,*),isblk1(4,*),id3x(*),i3x(*),                       11d14s23
     $     isblkk(4,*),jmats(*),kmats(*),itt(*),nbasdws(*)              10d30s23
      include "common.store"                                            7d21s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      igoul=itt(1)+3+nbasdws(1)*0
      zum0=0d0
      do is1=1,nsdlk1                                                   7d25s23
       run=0d0
       isb=isblk1(4,is1)                                                7d25s23
       if(isblk1(1,is1).eq.isblk1(2,is1))then                           7d25s23
        nrow=(nvirt(isblk1(1,is1))*(nvirt(isblk1(1,is1))+1))/2          11d14s23
        iswitch=0                                                       7d25s23
       else                                                             7d25s23
        nrow=nvirt(isblk1(1,is1))*nvirt(isblk1(2,is1))                  11d14s23
        iswitch=1                                                       7d25s23
       end if                                                           7d25s23
       call ilimts(irefo(isblk1(3,is1)),nvirt(isblk1(4,is1)),           11d14s23
     $      mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                   11d14s23
       nhere=ih+1-il                                                    11d15s23
       ncol=irefo(isblk1(3,is1))*nvirt(isb)                             7d25s23
       idd3=ibcoff                                                       11d14s23
       ibcoff=idd3+nrow*ncol                                             11d14s23
       call enough('dorbd3x.d3',bc,ibc)                                 11d14s23
       do iz=idd3,ibcoff-1                                               11d14s23
        bc(iz)=0d0                                                      11d14s23
       end do                                                           11d14s23
       jd3x=id3x(is1)                                                   11d14s23
       i10=i1s                                                          11d14s23
       i1n=irefo(isblk1(3,is1))                                         11d14s23
       do i2=i2s,i2e                                                    11d14s23
        if(i2.eq.i2e)i1n=i1e                                            11d14s23
        do i1=i10,i1n                                                   11d14s23
         jj=idd3+nrow*(i1-1+irefo(isblk1(3,is1))*(i2-1))                11d14s23
         do i34=0,nrow-1                                                11d14s23
          bc(jj+i34)=bc(jd3x+i34)                                       11d14s23
         end do                                                         11d14s23
         jd3x=jd3x+nrow                                                 11d14s23
        end do                                                          11d14s23
        i10=1                                                           11d14s23
       end do                                                           11d14s23
       call dws_gsumf(bc(idd3),nrow*ncol)                               11d14s23
c
       call ilimts(noc(isblk1(3,is1)),nvirt(isblk1(4,is1)),             11d15s23
     $      mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                   11d14s23
       nhere=ih+1-il                                                    11d15s23
       if(nhere.gt.0)then                                               11d15s23
        lsb=isblk1(1,is1)                                               11d15s23
        nmul=nhere*nvirt(isblk1(2,is1))                                 11d15s23
        itmpd=ibcoff                                                    11d15s23
        itmpi=itmpd+nvirt(lsb)*nmul                                     11d15s23
        ibcoff=itmpi+nmul*nvirt(lsb)                                    11d15s23
        call enough('dorbd3x.3xa',bc,ibc)                               11d15s23
        do iz=itmpd,ibcoff-1                                            11d15s23
         bc(iz)=0d0                                                     11d15s23
        end do                                                          11d15s23
        i10=i1s                                                         11d15s23
        i1n=noc(isblk1(3,is1))                                          11d15s23
        icol=0                                                          11d15s23
        nn=nvirt(isblk1(2,is1))*nhere                                   11d15s23
        ii=i3x(is1)                                                     11d15s23
        do i2=i2s,i2e                                                   11d15s23
         i2m=i2-1                                                       11d15s23
         if(i2.eq.i2e)i1n=i1e                                           11d15s23
         do i1=i10,i1n                                                  11d15s23
          if(i1.gt.idoubo(isblk1(3,is1)))then                           11d15s23
           i1m=i1-idoubo(isblk1(3,is1))-1                               11d15s23
           jdd=idd3+nrow*(i1m+irefo(isblk1(3,is1))*i2m)                 11d15s23
           do i4=0,nvirt(isblk1(2,is1))-1                                11d15s23
            jtmpd=itmpd+nvirt(lsb)*(i4+nvirt(isblk1(2,is1))*icol)        11d15s23
            jtmpi=itmpi+i4+nvirt(isblk1(2,is1))*icol                    11d15s23
            do i3=0,nvirt(lsb)-1                                         11d15s23
             irec=i3+nvirt(lsb)*i4                                       11d15s23
             ix=max(i3,i4)                                               11d15s23
             in=min(i3,i4)                                               11d15s23
             itri=((ix*(ix+1))/2)+in                                     11d15s23
             itri=itri+iswitch*(irec-itri)                               11d15s23
             bc(jtmpd+i3)=bc(itri+jdd)                                   11d15s23
             bc(jtmpi)=bc(itri+ii)                                       11d15s23
             jtmpi=jtmpi+nn                                              11d15s23
            end do                                                       11d15s23
            if(iswitch.eq.0)then                                         11d15s23
             itri=((i4*(i4+1))/2)+i4                                     11d15s23
             bc(jtmpd+i4)=bc(jtmpd+i4)+bc(itri+jdd)                      11d15s23
            end if                                                       11d15s23
           end do                                                        11d15s23
           icol=icol+1                                                   11d15s23
          end if                                                        11d15s23
          ii=ii+nrow                                                    11d15s23
         end do                                                         11d15s23
         i10=1                                                          11d15s23
        end do                                                          11d15s23
        jtt=itt(lsb)+noc(lsb)*(nbasdws(lsb)+1)                          11d15s23
        call dgemm('n','n',nvirt(lsb),nvirt(lsb),nmul,1d0,              11d15s23
     $       bc(itmpd),nvirt(lsb),bc(itmpi),nmul,1d0,                   11d15s23
     $       bc(jtt),nbasdws(lsb),'dorbd3x.3xa')                        11d15s23
        ibcoff=itmpd                                                    11d15s23
        if(isblk1(1,is1).ne.isblk1(2,is1))then                          11d15s23
         lsb=isblk1(2,is1)                                              11d15s23
         nmul=nhere*nvirt(isblk1(1,is1))                                11d15s23
         itmpd=ibcoff                                                    11d15s23
         itmpi=itmpd+nvirt(lsb)*nmul                                     11d15s23
         ibcoff=itmpi+nmul*nvirt(lsb)                                    11d15s23
         call enough('dorbd3x.3xb',bc,ibc)                               11d15s23
         do iz=itmpd,ibcoff-1                                            11d15s23
          bc(iz)=0d0                                                     11d15s23
         end do                                                          11d15s23
         i10=i1s                                                         11d15s23
         i1n=noc(isblk1(3,is1))                                          11d15s23
         icol=0                                                          11d15s23
         nn=nvirt(lsb)*nvirt(isblk1(1,is1))                             11d15s23
         ii=i3x(is1)                                                     11d15s23
         do i2=i2s,i2e                                                   11d15s23
          i2m=i2-1                                                       11d15s23
          if(i2.eq.i2e)i1n=i1e                                           11d15s23
          do i1=i10,i1n                                                  11d15s23
           if(i1.gt.idoubo(isblk1(3,is1)))then                           11d15s23
            i1m=i1-idoubo(isblk1(3,is1))-1                               11d15s23
            jdd=idd3+nrow*(i1m+irefo(isblk1(3,is1))*i2m)                 11d15s23
            do i4=0,nvirt(lsb)-1                                        11d15s23
             jtmpd=itmpd+i4+nn*icol                                     11d15s23
             jtmpi=itmpi+nvirt(isblk1(1,is1))*(icol+nhere*i4)           11d15s23
             do i3=0,nvirt(isblk1(1,is1))-1                             11d15s23
              irec=i3+nvirt(isblk1(1,is1))*i4                           11d15s23
              bc(jtmpd)=bc(irec+jdd)                                    11d15s23
              bc(jtmpi+i3)=bc(irec+ii)                                  11d15s23
              jtmpd=jtmpd+nvirt(lsb)
             end do                                                       11d15s23
            end do                                                        11d15s23
            icol=icol+1                                                   11d15s23
           end if                                                        11d15s23
           ii=ii+nrow                                                    11d15s23
          end do                                                         11d15s23
          i10=1                                                          11d15s23
         end do                                                          11d15s23
         jtt=itt(lsb)+noc(lsb)*(nbasdws(lsb)+1)                          11d15s23
         call dgemm('n','n',nvirt(lsb),nvirt(lsb),nmul,1d0,              11d15s23
     $       bc(itmpd),nvirt(lsb),bc(itmpi),nmul,1d0,                   11d15s23
     $       bc(jtt),nbasdws(lsb),'dorbd3x.3xb')                        11d15s23
         ibcoff=itmpd                                                   11d15s23
        end if                                                          11d15s23
        nv2=i2e+1-i2s                                                   11d15s23
        lsb=isblk1(3,is1)                                               11d15s23
        nmul=nrow*nv2                                                   11d15s23
        itmpd=ibcoff                                                    11d15s23
        itmpi=itmpd+irefo(lsb)*nmul                                     11d15s23
        ibcoff=itmpi+nmul*noc(lsb)                                      11d15s23
        call enough('dorbd3x.3xc',bc,ibc)                               11d15s23
        do iz=itmpd,ibcoff-1                                            11d15s23
         bc(iz)=0d0                                                     11d15s23
        end do                                                          11d15s23
        nn=irefo(lsb)*nrow                                              11d15s23
        do i2=i2s,i2e                                                   11d15s23
         i2m=i2-1                                                       11d15s23
         i2mm=i2-i2s                                                    11d15s23
         do i1=0,irefo(lsb)-1                                           11d15s23
          jdd3=idd3+nrow*(i1+irefo(lsb)*i2m)                            11d15s23
          jtmpd=itmpd+i1+nn*i2mm                                        11d15s23
          do i34=0,nrow-1                                               11d15s23
           bc(jtmpd)=bc(jdd3+i34)                                       11d15s23
           jtmpd=jtmpd+irefo(lsb)                                       11d15s23
          end do                                                        11d15s23
         end do                                                         11d15s23
        end do                                                          11d15s23
        i10=i1s                                                         11d15s23
        i1n=noc(lsb)                                                    11d15s23
        ii=i3x(is1)                                                     11d15s23
        do i2=i2s,i2e                                                   11d15s23
         i2m=i2-i2s                                                     11d15s23
         if(i2.eq.i2e)i1n=i1e                                           11d15s23
         do i1=i10,i1n                                                  11d15s23
          i1m=i1-1                                                      11d15s23
          jtmpi=itmpi+nrow*(i2m+nv2*i1m)                                11d15s23
          do i34=0,nrow-1                                               11d15s23
           bc(jtmpi+i34)=bc(ii+i34)                                     11d15s23
          end do                                                        11d15s23
          ii=ii+nrow                                                    11d15s23
         end do                                                         11d15s23
         i10=1                                                          11d15s23
        end do                                                          11d15s23
        jtt=itt(lsb)+idoubo(lsb)                                        11d15s23
        call dgemm('n','n',irefo(lsb),noc(lsb),nmul,1d0,                11d15s23
     $       bc(itmpd),irefo(lsb),bc(itmpi),nmul,1d0,                   11d15s23
     $       bc(jtt),nbasdws(lsb),'dorbd3x.3xc')                        11d15s23
        ibcoff=itmpd                                                    11d15s23
        lsb=isblk1(4,is1)                                               11d15s23
        nmul=nrow*irefo(isblk1(3,is1))                                  11d15s23
        itmpd=ibcoff                                                    11d15s23
        itmpi=itmpd+nvirt(lsb)*nmul                                     11d15s23
        ibcoff=itmpi+nmul*nv2                                           11d15s23
        call enough('dorbd3x.3xd',bc,ibc)                               11d15s23
        do iz=itmpd,ibcoff-1                                            11d15s23
         bc(iz)=0d0                                                     11d15s23
        end do                                                          11d15s23
        do i2=0,nvirt(lsb)-1                                            11d15s23
         do i1=0,nmul-1                                                 11d15s23
          i12=idd3+i1+nmul*i2                                           11d15s23
          i21=itmpd+i2+nvirt(lsb)*i1                                    11d15s23
          bc(i21)=bc(i12)                                               11d15s23
         end do                                                         11d15s23
        end do                                                          11d15s23
        i10=i1s                                                         11d15s23
        i1n=noc(isblk1(3,is1))                                          11d15s23
        ii=i3x(is1)                                                     11d15s23
        do i2=i2s,i2e                                                   11d15s23
         i2m=i2-i2s                                                     11d15s23
         if(i2.eq.i2e)i1n=i1e                                           11d15s23
         do i1=i10,i1n                                                  11d15s23
          if(i1.gt.idoubo(isblk1(3,is1)))then                           11d15s23
           i1m=i1-1-idoubo(isblk1(3,is1))                               11d15s23
           jtmpi=itmpi+nrow*(i1m+irefo(isblk1(3,is1))*i2m)              11d15s23
           do i34=0,nrow-1                                              11d15s23
            bc(jtmpi+i34)=bc(ii+i34)                                    11d15s23
           end do                                                       11d15s23
          end if                                                        11d15s23
          ii=ii+nrow                                                    11d15s23
         end do                                                         11d15s23
         i10=1                                                          11d15s23
        end do                                                          11d15s23
        jtt=itt(lsb)+noc(lsb)+nbasdws(lsb)*(noc(lsb)+i2s-1)             11d15s23
        call dgemm('n','n',nvirt(lsb),nv2,nmul,1d0,                     11d15s23
     $       bc(itmpd),nvirt(isb),bc(itmpi),nmul,1d0,                   11d15s23
     $       bc(jtt),nbasdws(lsb),'dorbd3x.3xd')                        11d15s23
        ibcoff=itmpd                                                    11d15s23
       end if                                                           11d15s23
       do is=1,nsdlk                                                    11d14s23
        if(isblk1(1,is1).eq.isblk(3,is).and.
     $     isblk1(2,is1).eq.isblk(4,is))then                            11d14s23
         if(isblk1(3,is1).eq.isblk(1,is))then                           11d14s23
          call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,   11d14s23
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          11d14s23
          if(iswitch.eq.0)then                                          11d14s23
           nrowj=(noc(isblk(1,is))*(noc(isblk(1,is))+1))/2              11d14s23
          else                                                          11d14s23
           nrowj=noc(isblk(1,is))*noc(isblk(2,is))                      11d14s23
          end if                                                        11d14s23
          nhere=ih+1-il                                                 11d14s23
          if(min(nhere,irefo(isblk(1,is)),noc(isb)).gt.0)then           11d14s23
           nmul=nhere*irefo(isblk(1,is))                                11d14s23
           itmpd=ibcoff                                                 11d14s23
           itmpi=itmpd+nvirt(isb)*nmul                                  11d14s23
           ibcoff=itmpi+nmul*noc(isb)                                   11d14s23
           call enough('dorbd3x.jd',bc,ibc)                             11d14s23
           do iz=itmpi,ibcoff-1                                         11d14s23
            bc(iz)=0d0                                                  11d14s23
           end do                                                       11d14s23
           i10=i1s                                                      11d14s23
           i1n=nvirt(isblk(3,is))                                       11d14s23
           icol=0                                                       11d14s23
           nn=nrow*irefo(isblk1(3,is1))                                 11d14s23
           nnn=nvirt(isb)*nhere                                         11d14s23
           nnnn=nvirt(isb)*irefo(isblk1(3,is1))                         11d14s23
           do i2=i2s,i2e                                                11d14s23
            i2m=i2-1                                                    11d14s23
            if(i2.eq.i2e)i1n=i1e                                        11d14s23
            do i1=i10,i1n                                               11d14s23
             i1m=i1-1                                                   11d14s23
             ix=max(i1m,i2m)                                            11d14s23
             in=min(i1m,i2m)                                            11d14s23
             irec=i1m+nvirt(isblk(3,is))*i2m                            11d14s23
             itri=((ix*(ix+1))/2)+in                                    11d14s23
             itri=itri+iswitch*(irec-itri)                              11d14s23
             ff=1d0                                                     11d14s23
             if(i1.ne.i2.and.iswitch.eq.0)ff=0.5d0                      11d14s23
             do i4=0,nvirt(isb)-1                                       11d14s23
              jdd3=idd3+itri+nn*i4                                      11d14s23
              jtmpd=itmpd+i4+nnnn*icol                                  11d14s23
              do i3=0,irefo(isblk1(3,is1))-1                            11d14s23
               bc(jtmpd)=bc(jdd3)*ff                                    11d14s23
               jdd3=jdd3+nrow                                           11d14s23
               jtmpd=jtmpd+nvirt(isb)                                   11d14s23
              end do                                                    11d14s23
             end do                                                     11d14s23
             icol=icol+1                                                11d14s23
            end do                                                      11d14s23
            i10=1                                                       11d14s23
           end do                                                       11d14s23
           i10=i1s                                                      11d14s23
           i1n=nvirt(isblk(3,is))                                       11d14s23
           jj=jmats(is)                                                 11d14s23
           icol=0                                                       11d14s23
           do i2=i2s,i2e                                                11d14s23
            if(i2.eq.i2e)i1n=i1e                                        11d14s23
            do i1=i10,i1n                                               11d14s23
             do i4=0,noc(isb)-1                                         11d14s23
              jtmpi=itmpi+irefo(isblk(1,is))*(icol+nhere*i4)            11d14s23
              do i3=0,irefo(isblk(1,is))-1                              11d14s23
               i3p=i3+idoubo(isblk(1,is))                               11d14s23
               irec=i3p+noc(isblk(1,is))*i4                             11d14s23
               ix=max(i3p,i4)                                           11d14s23
               in=min(i3p,i4)                                           11d14s23
               itri=((ix*(ix+1))/2)+in                                  11d14s23
               itri0=itri
               itri=itri+iswitch*(irec-itri)+jj                         11d14s23
               bc(jtmpi+i3)=bc(itri)                                    11d14s23
              end do                                                    11d14s23
             end do                                                     11d14s23
             jj=jj+nrowj                                                11d14s23
             icol=icol+1                                                11d14s23
            end do                                                      11d14s23
            i10=1                                                       11d14s23
           end do                                                       11d14s23
           jtt=itt(isb)+noc(isb)                                        11d14s23
           call dgemm('n','n',nvirt(isb),noc(isb),nmul,1d0,             11d14s23
     $          bc(itmpd),nvirt(isb),bc(itmpi),nmul,1d0,                11d14s23
     $          bc(jtt),nbasdws(isb),'dorbd3x.j')                       11d14s23
           ibcoff=itmpd                                                 11d14s23
          end if                                                        11d14s23
         else if(isblk1(3,is1).eq.isblk(2,is))then                      11d14s23
          call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,   11d14s23
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          11d14s23
          nrowj=noc(isblk(1,is))*noc(isblk(2,is))                       11d14s23
          nhere=ih+1-il                                                 11d14s23
          if(min(nhere,irefo(isblk(2,is)),noc(isb)).gt.0)then           11d14s23
           nmul=nhere*irefo(isblk(2,is))                                11d14s23
           itmpd=ibcoff                                                 11d14s23
           itmpi=itmpd+nvirt(isb)*nmul                                  11d14s23
           ibcoff=itmpi+nmul*noc(isb)                                   11d14s23
           call enough('dorbd3x.jdb',bc,ibc)                             11d14s23
           do iz=itmpi,ibcoff-1                                         11d14s23
            bc(iz)=0d0                                                  11d14s23
           end do                                                       11d14s23
           i10=i1s                                                      11d14s23
           i1n=nvirt(isblk(3,is))                                       11d14s23
           icol=0                                                       11d14s23
           nn=nrow*irefo(isblk1(3,is1))                                 11d14s23
           nnn=nvirt(isb)*nhere                                         11d14s23
           nnnn=nvirt(isb)*irefo(isblk1(3,is1))                         11d14s23
           do i2=i2s,i2e                                                11d14s23
            i2m=i2-1                                                    11d14s23
            if(i2.eq.i2e)i1n=i1e                                        11d14s23
            do i1=i10,i1n                                               11d14s23
             i1m=i1-1                                                   11d14s23
             irec=i1m+nvirt(isblk(3,is))*i2m                            11d14s23
             do i4=0,nvirt(isb)-1                                       11d14s23
              jdd3=idd3+irec+nn*i4                                      11d14s23
              jtmpd=itmpd+i4+nnnn*icol                                  11d14s23
              do i3=0,irefo(isblk1(3,is1))-1                            11d14s23
               bc(jtmpd)=bc(jdd3)                                       11d14s23
               jdd3=jdd3+nrow                                           11d14s23
               jtmpd=jtmpd+nvirt(isb)                                   11d14s23
              end do                                                    11d14s23
             end do                                                     11d14s23
             icol=icol+1                                                11d14s23
            end do                                                      11d14s23
            i10=1                                                       11d14s23
           end do                                                       11d14s23
           i10=i1s                                                      11d14s23
           i1n=nvirt(isblk(3,is))                                       11d14s23
           jj=jmats(is)                                                 11d14s23
           icol=0                                                       11d14s23
           do i2=i2s,i2e                                                11d14s23
            if(i2.eq.i2e)i1n=i1e                                        11d14s23
            do i1=i10,i1n                                               11d14s23
             do i4=0,irefo(isblk(2,is))-1                               11d14s23
              i4p=i4+idoubo(isblk(2,is))                                11d14s23
              jtmpi=itmpi+i4+irefo(isblk(2,is))*icol                    11d14s23
              do i3=0,noc(isb)-1                                        11d14s23
               irec=i3+noc(isblk(1,is))*i4p+jj                          11d14s23
               bc(jtmpi)=bc(irec)                                       11d14s23
               jtmpi=jtmpi+nmul                                         11d14s23
              end do                                                    11d14s23
             end do                                                     11d14s23
             jj=jj+nrowj                                                11d14s23
             icol=icol+1                                                11d14s23
            end do                                                      11d14s23
            i10=1                                                       11d14s23
           end do                                                       11d14s23
           jtt=itt(isb)+noc(isb)                                        11d14s23
           call dgemm('n','n',nvirt(isb),noc(isb),nmul,1d0,             11d14s23
     $          bc(itmpd),nvirt(isb),bc(itmpi),nmul,1d0,                11d14s23
     $          bc(jtt),nbasdws(isb),'dorbd3x.jb')                       11d14s23
           ibcoff=itmpd                                                 11d14s23
          end if                                                        11d14s23
         end if                                                         11d14s23
        end if                                                           11d14s23
       end do                                                           11d14s23
       do is=1,nsdlkk                                                   11d14s23
        if(isblkk(3,is).eq.isblk1(4,is1).and.                           11d14s23
     $       isblkk(2,is).eq.isblk1(3,is1))then                         11d14s23
         if(isblkk(4,is).eq.isblk1(2,is1).or.
     $      isblkk(4,is).eq.isblk1(1,is1))then                          11d14s23
          call ilimts(nvirt(isblkk(3,is)),nvirt(isblkk(4,is)),mynprocg, 11d14s23
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          11d14s23
          nhere=ih+1-il                                                 11d14s23
          if(min(nhere,irefo(isblk1(3,is1))).gt.0)then                  11d14s23
           nmul=nhere*irefo(isblk1(3,is1))                              11d14s23
           if(isblkk(4,is).eq.isblk1(2,is1))then                        11d14s23
            isb=isblk1(1,is1)                                           11d14s23
           else                                                         11d14s23
            isb=isblk1(2,is1)                                           11d14s23
           end if                                                       11d14s23
           if(noc(isb).gt.0)then                                        11d14s23
            itmpd=ibcoff                                                 11d14s23
            itmpi=itmpd+nvirt(isb)*nmul                                  11d14s23
            ibcoff=itmpi+nmul*noc(isb)                                   11d14s23
            call enough('dorbd3x.ka',bc,ibc)                             11d14s23
            do iz=itmpi,ibcoff-1                                         11d14s23
             bc(iz)=0d0                                                  11d14s23
            end do                                                       11d14s23
            ff=1d0                                                      11d14s23
            if(iswitch.eq.0)ff=0.5d0                                    11d14s23
            i10=i1s                                                      11d14s23
            i1n=nvirt(isblkk(3,is))                                      11d14s23
            icol=0                                                       11d14s23
            do i2=i2s,i2e                                                11d14s23
             i2m=i2-1                                                    11d14s23
             if(i2.eq.i2e)i1n=i1e                                        11d14s23
             if(isblkk(4,is).eq.isblk1(2,is1))then                       11d14s23
              do i1=i10,i1n                                              11d14s23
               i1m=i1-1                                                  11d14s23
               do i4=0,irefo(isblk1(3,is1))-1                            11d14s23
                jtmpd=itmpd+nvirt(isb)*(i4+irefo(isblk1(3,is1))*icol)    11d14s23
                jj=idd3+nrow*(i4+irefo(isblk1(3,is1))*i1m)               11d14s23
                do i3=0,nvirt(isb)-1                                     11d14s23
                 irec=i3+nvirt(isb)*i2m                                  11d14s23
                 ix=max(i3,i2m)                                          11d14s23
                 in=min(i3,i2m)                                          11d14s23
                 itri=((ix*(ix+1))/2)+in                                 11d14s23
                 itri=itri+iswitch*(irec-itri)+jj                        11d14s23
                 bc(jtmpd+i3)=bc(itri)                                  11d14s23
                end do                                                   11d14s23
                if(iswitch.eq.0)then                                     11d14s23
                 itri=((i2m*(i2m+1))/2)+i2m+jj                           11d14s23
                 bc(jtmpd+i2m)=bc(jtmpd+i2m)+bc(itri)                   11d14s23
                end if                                                   11d14s23
               end do                                                    11d14s23
               icol=icol+1                                               11d14s23
              end do                                                     11d14s23
             else                                                        11d14s23
              do i1=i10,i1n                                              11d14s23
               i1m=i1-1                                                  11d14s23
               do i4=0,irefo(isblk1(3,is1))-1                            11d14s23
                jtmpd=itmpd+nvirt(isb)*(i4+irefo(isblk1(3,is1))*icol)    11d14s23
                jj=idd3+nrow*(i4+irefo(isblk1(3,is1))*i1m)               11d14s23
                do i3=0,nvirt(isb)-1                                     11d14s23
                 irec=i2m+nvirt(isblkk(4,is))*i3+jj                      11d14s23
                 bc(jtmpd+i3)=bc(irec)                                   11d14s23
                end do                                                   11d14s23
               end do                                                    11d14s23
               icol=icol+1                                               11d14s23
              end do                                                     11d14s23
             end if                                                      11d14s23
             i10=1                                                       11d14s23
            end do                                                       11d14s23
            i10=i1s                                                      11d14s23
            i1n=nvirt(isblkk(3,is))                                      11d14s23
            kk=kmats(is)                                                 11d14s23
            nrowk=noc(isblkk(1,is))*noc(isblkk(2,is))                    11d14s23
            icol=0                                                       11d14s23
            do i2=i2s,i2e                                                11d14s23
             if(i2.eq.i2e)i1n=i1e                                        11d14s23
             do i1=i10,i1n                                               11d14s23
              do i4=0,irefo(isblkk(2,is))-1                              11d14s23
               i4p=i4+idoubo(isblkk(2,is))                               11d14s23
               kkk=kk+noc(isb)*i4p                                       11d14s23
               jtmpi=itmpi+i4+irefo(isblkk(2,is))*icol
               do i3=0,noc(isb)-1                                        11d14s23
                bc(jtmpi)=bc(kkk+i3)                                     11d14s23
                jtmpi=jtmpi+nmul                                         11d14s23
               end do                                                    11d14s23
              end do                                                     11d14s23
              icol=icol+1                                                11d14s23
              kk=kk+nrowk                                                11d14s23
             end do                                                      11d14s23
             i10=1                                                      11d15s23
            end do                                                       11d14s23
            jtt=itt(isb)+noc(isb)                                        11d14s23
            call dgemm('n','n',nvirt(isb),noc(isb),nmul,1d0,             11d14s23
     $           bc(itmpd),nvirt(isb),bc(itmpi),nmul,1d0,                11d14s23
     $          bc(jtt),nbasdws(isb),'dorbd3x.k')                       11d14s23
            ibcoff=itmpd                                                11d14s23
           end if                                                       11d14s23
          end if                                                        11d14s23
         end if                                                         11d14s23
        end if                                                          11d14s23
       end do                                                           11d14s23
       ibcoff=idd3                                                      11d14s23
      end do                                                            7d25s23
      return
      end
