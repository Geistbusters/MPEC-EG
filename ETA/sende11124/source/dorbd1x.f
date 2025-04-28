c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine dorbd1x(nsymb,idoubo,irefo,noc,nvirt,nsdlk,isblk,      6d24s24
     $     nsdlk1,isblk1,nsdlkk,isblkk,id1x,ioooo,ionex,jmats,kmats,    7d25s23
     $     bc,ibc,itt,nbasdws)                                          6d14s24
      implicit real*8 (a-h,o-z)
      dimension idoubo(*),irefo(*),noc(*),nvirt(*),                     6d24s24
     $     isblk(4,*),isblk1(4,*),id1x(*),ioooo(*),ionex(*),            6d14s24
     $     isblkk(4,*),jmats(*),kmats(*),itt(*),nbasdws(*)              10d30s23
      include "common.store"                                            7d21s23
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      igoul=itt(1)+1+nbasdws(1)*3
      zum0=0d0
      do is1=1,nsdlk1                                                   7d25s23
       isb=isblk1(4,is1)                                                7d25s23
       if(isblk1(1,is1).eq.isblk1(2,is1))then                           7d25s23
        nrow1=(irefo(isblk1(1,is1))*(irefo(isblk1(1,is1))+1))/2         7d25s23
        nrow2=(noc(isblk1(1,is1))*(noc(isblk1(1,is1))+1))/2             7d25s23
        iswitch=0                                                       7d25s23
       else                                                             7d25s23
        nrow1=irefo(isblk1(1,is1))*irefo(isblk1(2,is1))                 7d25s23
        nrow2=noc(isblk1(1,is1))*noc(isblk1(2,is1))                     7d25s23
        iswitch=1                                                       7d25s23
       end if                                                           7d25s23
       ncol=irefo(isblk1(3,is1))*nvirt(isb)                             7d25s23
       if(min(ncol,nrow1).gt.0)then                                     7d25s23
        nmul=nrow1*irefo(isblk1(3,is1))                                 7d25s23
        itmpd=ibcoff                                                    7d25s23
        ibcoff=itmpd+nmul*nvirt(isb)                                    7d25s23
        call enough('dorbd1x.tmpda',bc,ibc)                             7d25s23
        iad=id1x(is1)                                                   7d25s23
        do i3=0,irefo(isblk1(3,is1))-1                                  7d25s23
         do i12=0,nrow1-1                                               7d25s23
          jtmpd=itmpd+i12+nrow1*i3                                      7d25s23
          do i4=0,nvirt(isb)-1                                          7d25s23
           bc(jtmpd)=bc(iad+i4)                                         10d18s23
           jtmpd=jtmpd+nmul                                             7d25s23
          end do                                                        7d25s23
          iad=iad+nvirt(isb)                                            7d25s23
         end do                                                         7d25s23
        end do                                                          7d25s23
        do i=0,nmul*nvirt(isb)-1                                        7d25s23
         bc(id1x(is1)+i)=bc(itmpd+i)                                    7d25s23
        end do                                                          7d25s23
        ibcoff=itmpd                                                    7d25s23
        imap=ibcoff                                                     7d25s23
        ibcoff=imap+nrow1                                               7d25s23
        call enough('dorbd1x.imap',bc,ibc)
        jmap=imap                                                       7d25s23
        do i2=0,irefo(isblk1(2,is1))-1                                  7d25s23
         i2p=i2+idoubo(isblk1(2,is1))                                   7d25s23
         itop=i2+iswitch*(irefo(isblk1(1,is1))-1-i2)                    7d25s23
         irec=idoubo(isblk1(1,is1))+noc(isblk1(1,is1))*i2p                               7d25s23
         itri=((i2p*(i2p+1))/2)+idoubo(isblk1(1,is1))                   7d25s23
         itri=itri+iswitch*(irec-itri)                                  7d25s23
         do i1=0,itop                                                   7d25s23
          ibc(jmap)=itri+i1                                             7d25s23
          jmap=jmap+1                                                   7d25s23
         end do                                                         7d25s23
        end do                                                          7d25s23
        jsb=isblk1(3,is1)                                               7d25s23
        do is=1,nsdlk                                                   7d25s23
         if(isblk1(1,is1).eq.isblk(1,is).and.                           7d25s23
     $      isblk1(2,is1).eq.isblk(2,is).and.                           7d25s23
     $      isblk1(3,is1).eq.isblk(3,is))then                           7d25s23
          call ilimts(nvirt(isblk(3,is)),nvirt(isblk(4,is)),mynprocg,   7d27s23
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          7d27s23
          nhere=ih+1-il                                                 7d27s23
          if(nhere.gt.0)then                                            7d27s23
           nv=i2e+1-i2s                                                 7d27s23
           if(irefo(jsb).gt.0)then                                      7d28s23
            nmult=nrow1*nv                                              7d27s23
            itmpd=ibcoff                                                7d27s23
            itmpi=itmpd+irefo(jsb)*nmult                                7d28s23
            ibcoff=itmpi+nmult*nvirt(jsb)                               7d27s23
            call enough('dorbd1x.ev',bc,ibc)                            7d27s23
            do iz=itmpi,ibcoff-1                                        7d27s23
             bc(iz)=0d0                                                 7d27s23
            end do                                                      7d27s23
            do iv=i2s,i2e                                               7d27s23
             ivm=iv-1                                                   7d27s23
             ivmm=iv-i2s                                                7d27s23
             do i3=0,irefo(jsb)-1                                       7d28s23
              iadd=id1x(is1)+nrow1*(i3+irefo(jsb)*ivm)                  7d27s23
              iad=itmpd+i3+irefo(jsb)*nrow1*ivmm                        7d28s23
              do i12=0,nrow1-1                                          7d27s23
               bc(iad)=bc(iadd+i12)                                     7d27s23
               iad=iad+irefo(jsb)                                       7d28s23
              end do                                                    7d27s23
             end do                                                     7d27s23
            end do                                                      7d27s23
            i10=i1s                                                     7d27s23
            i1n=nvirt(jsb)                                              7d27s23
            ii=jmats(is)                                                7d27s23
            do i2=i2s,i2e                                               7d27s23
             i2m=i2-i2s                                                 7d27s23
             if(i2.eq.i2e)i1n=i1e                                       7d27s23
             do i1=i10,i1n                                              7d27s23
              i1m=i1-1                                                  7d27s23
              jtmpi=itmpi+nrow1*(i2m+nv*i1m)                            7d27s23
              do i34=0,nrow1-1                                          7d27s23
               iii=ii+ibc(imap+i34)                                     7d27s23
               bc(jtmpi+i34)=bc(iii)                                    7d27s23
              end do                                                    7d27s23
              ii=ii+nrow2                                               7d27s23
             end do                                                     7d27s23
             i10=1                                                      7d27s23
            end do                                                      7d27s23
            jtt=itt(jsb)+idoubo(jsb)+nbasdws(jsb)*noc(jsb)              10d30s23
            call dgemm('n','n',irefo(jsb),nvirt(jsb),nmult,1d0,         7d28s23
     $           bc(itmpd),irefo(jsb),bc(itmpi),nmult,1d0,              10d30s23
     $           bc(jtt),nbasdws(jsb),'dorbd1x.ev')                     10d30s23
            ibcoff=itmpd                                                7d27s23
           end if                                                       7d27s23
          end if                                                        7d27s23
          call ilimts(noc(isblk(3,is)),noc(isblk(4,is)),mynprocg,       7d25s23
     $        mynowprog,il,ih,i1s,i1e,i2s,i2e)                          7d25s23
          nhere=ih+1-il                                                 7d25s23
          if(nhere.gt.0)then                                            7d25s23
           if(noc(isb).gt.0)then                                        10d30s23
            itmpi=ibcoff                                                  7d25s23
            ibcoff=itmpi+nmul*noc(isb)                                  10d30s23
            call enough('dorbd1x.tmpi4oa',bc,ibc)                        7d25s23
            do iz=itmpi,ibcoff-1                                         7d25s23
             bc(iz)=0d0                                                  7d25s23
            end do                                                       7d25s23
            jtmpi=itmpi                                                  7d25s23
            do idd=0,noc(isb)-1                                         10d30s23
             do i1=0,irefo(isblk1(3,is1))-1                              7d25s23
              i1p=i1+idoubo(isblk1(3,is1))                               7d25s23
              icol=i1p+1+noc(isblk1(3,is1))*idd                          7d25s23
              if(icol.ge.il.and.icol.le.ih)then                          7d25s23
               icol0=icol                                                7d25s23
               icol=ioooo(is)+nrow2*(icol-il)                            7d25s23
               iad=itmpi+idd+noc(isb)*nrow1*i1                          10d30s23
               do i=0,nrow1-1                                            7d25s23
                ii=icol+ibc(imap+i)                                      7d25s23
                bc(iad)=bc(ii)                                           7d25s23
                iad=iad+noc(isb)                                        10d30s23
               end do                                                    7d25s23
              end if                                                     7d25s23
             end do                                                      7d25s23
            end do                                                       7d25s23
            itmp2=ibcoff                                                10d30s23
            ibcoff=itmp2+noc(isb)*nvirt(isb)                            10d30s23
            call enough('dorbd1x.tmp2',bc,ibc)                          10d30s23
            call dgemm('n','n',noc(isb),nvirt(isb),nmul,1d0,            10d30s23
     $          bc(itmpi),noc(isb),bc(id1x(is1)),nmul,0d0,              10d30s23
     $          bc(itmp2),noc(isb),'dorbd1x.4oa')                       10d30s23
            jtmp2=itmp2                                                 10d30s23
            do iv=0,nvirt(isb)-1                                        10d30s23
             jtt=itt(isb)+noc(isb)+iv                                   10d30s23
             do i=0,noc(isb)-1                                          10d30s23
              bc(jtt)=bc(jtt)+bc(jtmp2+i)                               10d30s23
              jtt=jtt+nbasdws(isb)                                      10d30s23
             end do                                                     10d30s23
             jtmp2=jtmp2+noc(isb)                                       10d30s23
            end do                                                      10d30s23
            ibcoff=itmpi                                                10d30s23
           end if
          end if                                                        7d25s23
         end if                                                         7d25s23
        end do                                                          7d25s23
        call ilimts(noc(jsb),nvirt(isb),mynprocg,mynowprog,il,ih,i1s,   7d25s23
     $       i1e,i2s,i2e)                                               7d25s23
        nhere=ih+1-il                                                   10d30s23
        nv=i2e+1-i2s                                                    7d25s23
        ksb=isblk1(1,is1)                                               7d25s23
        lsb=isblk1(2,is1)                                               7d27s23
        if(min(nv,irefo(jsb),nrow1).gt.0)then                           10d30s23
         nmul=nrow1*irefo(jsb)                                          10d30s23
         itmpd=ibcoff                                                   10d30s23
         itmpi=itmpd+nvirt(isb)*nmul                                    10d30s23
         ibcoff=itmpi+nv*nmul                                           10d30s23
         call enough('do4bd1x.tmpdvv',bc,ibc)                           10d30s23
         do iz=itmpi,ibcoff-1                                           10d30s23
          bc(iz)=0d0                                                    10d30s23
         end do                                                         10d30s23
         do iv=0,nvirt(isb)-1                                           10d30s23
          do i123=0,nmul-1                                              10d31s23
           ij=id1x(is1)+i123+nmul*iv                                    10d30s23
           ji=itmpd+iv+nvirt(isb)*i123                                  10d31s23
           bc(ji)=bc(ij)                                                10d30s23
          end do                                                        10d30s23
         end do                                                         10d30s23
         ii=ionex(is1)                                                  10d30s23
         i10=i1s                                                        10d30s23
         i1n=noc(jsb)                                                   10d30s23
         do i2=i2s,i2e                                                  10d30s23
          if(i2.eq.i2e)i1n=i1e                                          10d30s23
          do i1=i10,i1n                                                 10d30s23
           if(i1.gt.idoubo(jsb))then                                    10d30s23
            jtmpi=itmpi+nrow1*(i1-1-idoubo(jsb)+irefo(jsb)*(i2-i2s))    10d30s23
            do i=0,nrow1-1                                               10d30s23
             iii=ii+ibc(imap+i)                                          10d30s23
             bc(jtmpi+i)=bc(iii)                                        10d31s23
            end do                                                       10d30s23
           end if                                                       10d30s23
           ii=ii+nrow2                                                  10d30s23
          end do                                                        10d30s23
          i10=1                                                         10d30s23
         end do                                                         10d30s23
         jtt=itt(isb)+noc(isb)+nbasdws(isb)*(noc(isb)+i2s-1)            10d30s23
         call dgemm('n','n',nvirt(isb),nv,nmul,1d0,                     10d30s23
     $        bc(itmpd),nvirt(isb),bc(itmpi),nmul,1d0,                  10d30s23
     $        bc(jtt),nbasdws(isb),'dorbd1x.goul1c')                    10d30s23
         ibcoff=itmpd                                                   10d30s23
        end if                                                          10d30s23
        if(min(nv,noc(ksb),irefo(lsb),irefo(jsb)).gt.0)then             10d30s23
         nmul=irefo(lsb)*irefo(jsb)*nv                                  7d27s23
         itmpd=ibcoff                                                   7d25s23
         itmpi=itmpd+irefo(ksb)*nmul                                    10d30s23
         ibcoff=itmpi+noc(ksb)*nmul                                     10d30s23
         call enough('dorbd1x.ea12',bc,ibc)                             7d25s23
         do iz=itmpi,ibcoff-1                                           7d25s23
          bc(iz)=0d0                                                    7d25s23
         end do                                                         7d25s23
         do i2=0,irefo(lsb)-1                                           7d27s23
          do i1=0,irefo(ksb)-1                                          10d30s23
           ff=1d0                                                       7d27s23
           if(iswitch.eq.0)then                                         7d27s23
            ix=max(i1,i2)                                               7d27s23
            in=min(i1,i2)                                               7d27s23
            itri=((ix*(ix+1))/2)+in                                     7d27s23
            if(i1.ne.i2)ff=0.5d0                                        7d27s23
           else                                                         7d27s23
            itri=i1+irefo(ksb)*i2                                       7d27s23
           end if                                                       7d27s23
           do iv=i2s,i2e                                                7d27s23
            ivm=iv-1                                                    7d27s23
            ivmm=iv-i2s                                                 7d27s23
            do i3=0,irefo(jsb)-1                                        7d27s23
             iadd=id1x(is1)+itri+nrow1*(i3+irefo(jsb)*ivm)              7d27s23
             iad=itmpd+i1+irefo(ksb)*(i2+irefo(lsb)*(i3                 10d30s23
     $            +irefo(jsb)*ivmm))                                    7d25s23
             bc(iad)=bc(iadd)*ff                                        7d25s23
            end do                                                      7d25s23
           end do                                                       7d25s23
          end do                                                        7d25s23
         end do                                                         7d25s23
         do iv=i2s,i2e                                                  7d25s23
          ivm=iv-1                                                      7d25s23
          ivmm=iv-i2s                                                   7d25s23
          do i1=0,irefo(jsb)-1                                          7d25s23
           i1p=i1+idoubo(jsb)                                           7d25s23
           icol=i1p+1+noc(jsb)*ivm                                      7d25s23
           if(icol.ge.il.and.icol.le.ih)then                            7d25s23
            icol=ionex(is1)+nrow2*(icol-il)                             7d25s23
            do i4=0,irefo(lsb)-1                                        7d27s23
             i4p=i4+idoubo(lsb)                                         7d27s23
             do i3=0,noc(ksb)-1                                         10d30s23
              irec=i3+noc(ksb)*i4p                                      10d30s23
              ix=max(i4p,i3)                                            10d30s23
              in=min(i4p,i3)                                            10d30s23
              itri=((ix*(ix+1))/2)+in                                   7d25s23
              itri=itri+iswitch*(irec-itri)+icol                        7d25s23
              jtmpi=itmpi+i4+irefo(lsb)*(i1                             7d27s23
     $             +irefo(jsb)*(ivmm+nv*i3))                            7d27s23
              bc(jtmpi)=bc(itri)                                        7d25s23
             end do                                                     7d25s23
            end do                                                      7d25s23
           end if                                                       7d25s23
          end do                                                        7d25s23
         end do                                                         7d25s23
         jtt=itt(ksb)+idoubo(ksb)                                       10d30s23
         if(isb.eq.jsb)then                                             10d31s23
          fmul=2d0                                                      10d31s23
         else                                                           10d31s23
          fmul=1d0                                                      10d31s23
         end if                                                         10d31s23
         call dgemm('n','n',irefo(ksb),noc(ksb),nmul,fmul,              10d31s23
     $        bc(itmpd),irefo(ksb),bc(itmpi),nmul,1d0,                  10d30s23
     $        bc(jtt),nbasdws(ksb),'dorbd1x.1xea')                      7d27s23
         ibcoff=itmpi                                                   10d30s23
        end if                                                          10d30s23
        if(ksb.ne.lsb.and.                                              10d30s23
     $       min(nv,noc(lsb),irefo(ksb),irefo(jsb)).gt.0)then           10d30s23
         nmul=irefo(ksb)*irefo(jsb)*nv                                  7d27s23
         itmpd=ibcoff                                                   7d25s23
         itmpi=itmpd+irefo(lsb)*nmul                                    10d30s23
         ibcoff=itmpi+noc(lsb)*nmul                                     10d30s23
         call enough('dorbd1x.ea12',bc,ibc)                             7d25s23
         do iz=itmpi,ibcoff-1                                           7d25s23
          bc(iz)=0d0                                                    7d25s23
         end do                                                         7d25s23
         do i2=0,irefo(lsb)-1                                           7d27s23
          do i1=0,irefo(ksb)-1                                          10d30s23
           irec=i1+irefo(ksb)*i2                                        10d31s23
           do iv=i2s,i2e                                                7d27s23
            ivm=iv-1                                                    7d27s23
            ivmm=iv-i2s                                                 7d27s23
            do i3=0,irefo(jsb)-1                                        7d27s23
             iadd=id1x(is1)+irec+nrow1*(i3+irefo(jsb)*ivm)              10d31s23
             iad=itmpd+i2+irefo(lsb)*(i1+irefo(ksb)*(i3                 10d31s23
     $            +irefo(jsb)*ivmm))                                    10d31s23
             bc(iad)=bc(iadd)                                           10d31s23
            end do                                                      7d25s23
           end do                                                       7d25s23
          end do                                                        7d25s23
         end do                                                         7d25s23
         do iv=i2s,i2e                                                  7d25s23
          ivm=iv-1                                                      7d25s23
          ivmm=iv-i2s                                                   7d25s23
          do i1=0,irefo(jsb)-1                                          7d25s23
           i1p=i1+idoubo(jsb)                                           7d25s23
           icol=i1p+1+noc(jsb)*ivm                                      7d25s23
           if(icol.ge.il.and.icol.le.ih)then                            7d25s23
            icol=ionex(is1)+nrow2*(icol-il)                             7d25s23
            do i4=0,noc(lsb)-1                                          10d31s23
             do i3=0,irefo(ksb)-1                                       10d31s23
              i3p=i3+idoubo(ksb)                                        10d31s23
              irec=i3p+noc(ksb)*i4+icol                                 10d31s23
              jtmpi=itmpi+i3+irefo(ksb)*(i1                             10d30s23
     $             +irefo(jsb)*(ivmm+nv*i4))                            10d30s23
              bc(jtmpi)=bc(irec)                                        10d31s23
             end do                                                     7d25s23
            end do                                                      7d25s23
           end if                                                       7d25s23
          end do                                                        7d25s23
         end do                                                         7d25s23
         jtt=itt(lsb)+idoubo(lsb)                                       10d30s23
         call dgemm('n','n',irefo(lsb),noc(lsb),nmul,1d0,               10d31s23
     $        bc(itmpd),irefo(lsb),bc(itmpi),nmul,1d0,                  10d30s23
     $        bc(jtt),nbasdws(lsb),'dorbd1x.1xea')                      10d30s23
         ibcoff=itmpi                                                   10d30s23
        end if                                                          10d30s23
        if(min(nv,noc(jsb),irefo(jsb)).gt.0)then                        10d30s23
         nmul=nrow1*nv                                                  7d25s23
         itmpd=ibcoff                                                   7d25s23
         itmpi=itmpd+nmul*irefo(jsb)                                    10d30s23
         ibcoff=itmpi+nmul*noc(jsb)                                     10d30s23
         call enough('dorbd1x.1xda',bc,ibc)                             7d25s23
         do iz=itmpi,ibcoff-1                                           7d25s23
          bc(iz)=0d0                                                    7d25s23
         end do                                                         7d25s23
         do iv=i2s,i2e                                                  7d25s23
          ivm=iv-1                                                      7d25s23
          ivmm=iv-i2s                                                   7d25s23
          do ia=0,irefo(jsb)-1                                          10d30s23
           iadd=id1x(is1)+nrow1*(ia+irefo(jsb)*ivm)                     10d31s23
           jtmpd=itmpd+nrow1*(ivmm+nv*ia)                               7d25s23
           do i12=0,nrow1-1                                             7d25s23
            bc(jtmpd+i12)=bc(iadd+i12)                                  7d25s23
           end do                                                       7d25s23
          end do                                                        7d25s23
         end do                                                         7d25s23
         do iv=i2s,i2e                                                  7d25s23
          ivm=iv-1                                                      7d25s23
          ivmm=iv-i2s                                                   7d25s23
          do idd=0,noc(jsb)-1                                           10d30s23
           icol=idd+1+noc(jsb)*ivm                                       7d25s23
           if(icol.ge.il.and.icol.le.ih)then                            7d25s23
            icol=ionex(is1)+nrow2*(icol-il)                             7d25s23
            jtmpi=itmpi+idd+noc(jsb)*nrow1*ivmm                         10d30s23
            do i34=0,nrow1-1                                            7d25s23
             ii=icol+ibc(imap+i34)                                      7d25s23
             bc(jtmpi)=bc(ii)                                           7d25s23
             jtmpi=jtmpi+noc(jsb)                                       10d30s23
            end do                                                      7d25s23
           end if                                                       7d25s23
          end do                                                        7d25s23
         end do                                                         7d25s23
         itmp2=ibcoff                                                   10d30s23
         ibcoff=itmp2+noc(jsb)*irefo(jsb)                               10d30s23
         call enough('dorbd1x.tmp21b',bc,ibc)                           10d30s23
         call dgemm('n','n',noc(jsb),irefo(jsb),nmul,1d0,               10d30s23
     $        bc(itmpi),noc(jsb),bc(itmpd),nmul,0d0,                    10d30s23
     $        bc(itmp2),noc(jsb),'dorbd1x.da')                          10d30s23
         jtmp2=itmp2                                                    10d30s23
         do ia=0,irefo(jsb)-1                                           10d30s23
          jtt=itt(jsb)+idoubo(jsb)+ia                                   10d30s23
          do i=0,noc(jsb)-1                                             10d30s23
           bc(jtt)=bc(jtt)+bc(jtmp2+i)                                  10d30s23
           jtt=jtt+nbasdws(jsb)                                         10d30s23
          end do                                                        10d30s23
          jtmp2=jtmp2+noc(jsb)                                          10d30s23
         end do                                                         10d30s23
         ibcoff=itmpi                                                   10d30s23
        end if                                                          10d30s23
        isb=isblk1(1,is1)                                               7d28s23
        lsb=isblk1(2,is1)                                               7d31s23
        if(irefo(isb).gt.0)then                                         7d28s23
         do isk=1,nsdlkk                                                 7d28s23
          if(isblk1(2,is1).eq.isblkk(1,isk).and.                         7d28s23
     $      isblk1(3,is1).eq.isblkk(2,isk).and.                         7d28s23
     $      isblk1(4,is1).eq.isblkk(3,isk))then                         7d28s23
           call ilimts(nvirt(isblkk(3,isk)),nvirt(isblkk(4,isk)),        7d28s23
     $        mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                 7d28s23
           nhere=ih+1-il                                                 7d28s23
           if(nhere.gt.0)then                                            7d28s23
            nrowk=noc(isblkk(1,isk))*noc(isblkk(2,isk))                  7d28s23
            nrowi=irefo(isblkk(1,isk))*irefo(isblkk(2,isk))              7d28s23
            nv=i2e+1-i2s                                                 7d28s23
            nmult=nrowi*nvirt(isblkk(3,isk))                            7d31s23
            imapk=ibcoff                                                 7d28s23
            itmpd=imapk+nrowi                                            7d28s23
            itmpi=itmpd+nmult*irefo(isb)                                 7d28s23
            ibcoff=itmpi+nmult*nv                                       7d31s23
            call enough('dorbd1x.ktmpd',bc,ibc)                          7d28s23
            do iz=itmpi,ibcoff-1                                         7d28s23
             bc(iz)=0d0                                                  7d28s23
            end do                                                       7d28s23
            jmapk=imapk                                                  7d28s23
            do i2=0,irefo(isblkk(2,isk))-1                               7d28s23
             i2p=i2+idoubo(isblkk(2,isk))                                7d28s23
             do i1=0,irefo(isblkk(1,isk))-1                              7d28s23
              i1p=i1+idoubo(isblkk(1,isk))                               7d28s23
              ibc(jmapk+i1)=i1p+noc(isblkk(1,isk))*i2p                   7d28s23
             end do                                                      7d28s23
             jmapk=jmapk+irefo(isblkk(1,isk))                            7d28s23
            end do                                                       7d28s23
            nn=irefo(isblk1(1,is1))*irefo(isblk1(2,is1))                 7d28s23
            do i2=0,irefo(isblk1(2,is1))-1                               7d28s23
             do i1=0,irefo(isblk1(1,is1))-1                              7d28s23
              irec=i1+irefo(isblk1(1,is1))*i2                            7d28s23
              ix=max(i1,i2)                                              7d28s23
              in=min(i1,i2)                                              7d28s23
              itri=((ix*(ix+1))/2)+in                                    7d28s23
              itri=itri+iswitch*(irec-itri)                              7d28s23
              ff=0.5d0                                                   7d28s23
              if(iswitch.eq.0.and.i1.eq.i2)ff=1d0                        7d28s23
              do iv=0,nvirt(isblkk(3,isk))-1                            7d31s23
               iadd=id1x(is1)+itri+nrow1*irefo(isblk1(3,is1))*iv        7d31s23
               jtmpd=itmpd+i1+irefo(isblk1(1,is1))*(i2
     $              +irefo(isblk1(2,is1))*irefo(isblk1(3,is1))*iv)      7d31s23
               do i3=0,irefo(isblk1(3,is1))-1                            7d28s23
                bc(jtmpd)=bc(iadd)*ff                                    7d28s23
                jtmpd=jtmpd+nn                                           7d28s23
                iadd=iadd+nrow1                                          7d28s23
               end do                                                    7d28s23
              end do                                                     7d28s23
             end do                                                      7d28s23
            end do                                                       7d28s23
            i10=i1s                                                      7d28s23
            i1n=nvirt(isblkk(3,isk))                                     7d28s23
            kk=kmats(isk)                                                7d28s23
            do i2=i2s,i2e                                                7d28s23
             if(i2.eq.i2e)i1n=i1e                                        7d28s23
             do i1=i10,i1n                                               7d28s23
              jtmpi=itmpi+nrowi*(i1-1+nvirt(isblkk(3,isk))*(i2-i2s))    7d31s23
              do i34=0,nrowi-1                                           7d28s23
               ii=kk+ibc(imapk+i34)                                      7d28s23
               bc(jtmpi+i34)=bc(ii)                                      7d28s23
              end do                                                     7d28s23
              kk=kk+nrowk                                                7d28s23
             end do                                                      7d28s23
             i10=1                                                       7d28s23
            end do                                                       7d28s23
            jtt=itt(isb)+idoubo(isb)+nbasdws(isb)*(noc(isb)+i2s-1)      5d17s24
            call dgemm('n','n',irefo(isb),nv,nmult,2d0,                 7d31s23
     $          bc(itmpd),irefo(isb),bc(itmpi),nmult,1d0,               10d30s23
     $          bc(jtt),nbasdws(isb),'dorbd1x,tmpk')                     7d28s23
           end if                                                        7d28s23
          else if(isblk1(1,is1).ne.isblk1(2,is1).and.isblk1(1,is1).eq.
     $          isblkk(1,isk).and.                                      7d31s23
     $      isblk1(3,is1).eq.isblkk(2,isk).and.                         7d28s23
     $      isblk1(4,is1).eq.isblkk(3,isk))then                         7d28s23
           call ilimts(nvirt(isblkk(3,isk)),nvirt(isblkk(4,isk)),        7d28s23
     $        mynprocg,mynowprog,il,ih,i1s,i1e,i2s,i2e)                 7d28s23
           nhere=ih+1-il                                                 7d28s23
           if(nhere.gt.0)then                                            7d28s23
            nrowk=noc(isblkk(1,isk))*noc(isblkk(2,isk))                  7d28s23
            nrowi=irefo(isblkk(1,isk))*irefo(isblkk(2,isk))              7d28s23
            nv=i2e+1-i2s                                                 7d28s23
            nmult=nrowi*nvirt(isblkk(3,isk))                            7d31s23
            imapk=ibcoff                                                 7d28s23
            itmpd=imapk+nrowi                                            7d28s23
            itmpi=itmpd+nmult*irefo(lsb)                                 7d28s23
            ibcoff=itmpi+nmult*nv                                       7d31s23
            call enough('dorbd1x.ktmpdb',bc,ibc)                          7d28s23
            do iz=itmpi,ibcoff-1                                         7d28s23
             bc(iz)=0d0                                                  7d28s23
            end do                                                       7d28s23
            jmapk=imapk                                                  7d28s23
            do i2=0,irefo(isblkk(2,isk))-1                               7d28s23
             i2p=i2+idoubo(isblkk(2,isk))                                7d28s23
             do i1=0,irefo(isblkk(1,isk))-1                              7d28s23
              i1p=i1+idoubo(isblkk(1,isk))                               7d28s23
              ibc(jmapk+i1)=i1p+noc(isblkk(1,isk))*i2p                   7d28s23
             end do                                                      7d28s23
             jmapk=jmapk+irefo(isblkk(1,isk))                            7d28s23
            end do                                                       7d28s23
            nn=irefo(isblk1(1,is1))*irefo(isblk1(2,is1))                 7d28s23
            do i2=0,irefo(isblk1(2,is1))-1                               7d28s23
             do i1=0,irefo(isblk1(1,is1))-1                              7d28s23
              irec=i1+irefo(isblk1(1,is1))*i2                            7d28s23
              ff=0.5d0                                                   7d28s23
              do iv=0,nvirt(isblkk(3,isk))-1                            7d31s23
               iadd=id1x(is1)+irec+nrow1*irefo(isblk1(3,is1))*iv        7d31s23
               jtmpd=itmpd+i2+irefo(isblk1(2,is1))*(i1                  7d31s23
     $              +irefo(isblk1(1,is1))*irefo(isblk1(3,is1))*iv)      7d31s23
               do i3=0,irefo(isblk1(3,is1))-1                            7d28s23
                bc(jtmpd)=bc(iadd)*ff                                    7d28s23
                jtmpd=jtmpd+nn                                           7d28s23
                iadd=iadd+nrow1                                          7d28s23
               end do                                                    7d28s23
              end do                                                     7d28s23
             end do                                                      7d28s23
            end do                                                       7d28s23
            i10=i1s                                                      7d28s23
            i1n=nvirt(isblkk(3,isk))                                     7d28s23
            kk=kmats(isk)                                                7d28s23
            do i2=i2s,i2e                                                7d28s23
             if(i2.eq.i2e)i1n=i1e                                        7d28s23
             do i1=i10,i1n                                               7d28s23
              jtmpi=itmpi+nrowi*(i1-1+nvirt(isblkk(3,isk))*(i2-i2s))    7d31s23
              do i34=0,nrowi-1                                           7d28s23
               ii=kk+ibc(imapk+i34)                                      7d28s23
               bc(jtmpi+i34)=bc(ii)                                      7d28s23
              end do                                                     7d28s23
              kk=kk+nrowk                                                7d28s23
             end do                                                      7d28s23
             i10=1                                                       7d28s23
            end do                                                       7d28s23
            jtt=itt(lsb)+idoubo(lsb)+nbasdws(lsb)*(noc(lsb)+i2s-1)      5d20s24
            call dgemm('n','n',irefo(lsb),nv,nmult,2d0,                 7d31s23
     $          bc(itmpd),irefo(lsb),bc(itmpi),nmult,1d0,               10d30s23
     $          bc(jtt),nbasdws(lsb),'dorbd1x,tmpk')                    10d30s23
            ibcoff=imapk                                                10d31s23
           end if                                                       7d31s23
          end if                                                         7d28s23
         end do                                                          7d28s23
 7575    continue                                                        7d28s23
        end if                                                          7d28s23
        ibcoff=imap                                                     7d25s23
       end if                                                           7d25s23
      end do                                                            7d25s23
      return
      end
