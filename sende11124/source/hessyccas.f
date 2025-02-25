c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hessyccas(itrial,imt,iamat,ihessa,ihessb,ihessc,idoub, 12d23s17
     $     iacto,noc,nvirtc,nsymb,xlami,bc,ibc)                         11d10s22
      implicit real*8 (a-h,o-z)
c
c     form matrix vector product for augmented hessian optimization.    11d29s17
c     recall ibmat is half the gradiant.                                11d29s17
c     for now, hessa is the same for each proc, and ihessc is not       11d29s17
c     symmetrized. the next step has hessa not yet summed over procs and11d29s17
c     two parts of hessb. Then the only global summing in buildhesscas  11d29s17
c     will be for the dd, aa, and da parts.
c     we can also explicitly treat the dd, aa and da parts and thus     11d29s17
c     avoid all global sums in buildhesscas.                            11d29s17
c     at any rate, the plan is to compute this procs contribution to the11d29s17
c     matrix vector product, and leave all global sums to the calling   11d29s17
c     programs discression.                                             11d29s17
c
      include "common.store"
      dimension itrial(nsymb),imt(nsymb),ihessa(8,8),ihessb(8,8),
     $     ihessc(8,8),idoub(nsymb),iacto(nsymb),noc(nsymb),
     $     nvirtc(nsymb),iamat(nsymb)                                   12d23s17
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      save                                                              4d3s18
      data icall/0/                                                     4d3s18
c
c     the order of variables is active,double, then virt,occ
c
      const=bc(itrial(1))                                               11d29s17
      do isb=1,nsymb                                                    12d24s17
       nv=idoub(isb)*iacto(isb)+noc(isb)*nvirtc(isb)                    11d29s17
       if(isb.eq.1)nv=nv+1                                              11d29s17
       do i=0,nv-1                                                      11d29s17
        bc(imt(isb)+i)=0d0                                              11d29s17
       end do                                                           11d29s17
      end do                                                            12d24s17
      do isb=1,nsymb
       nv=idoub(isb)*iacto(isb)+noc(isb)*nvirtc(isb)                    11d29s17
       if(isb.eq.1)nv=nv+1                                              11d29s17
       sumb=0d0                                                         12d23s17
       do isk=1,nsymb                                                   11d29s17
c
c      gradient dot trial vector which goes into constraint
c
        if(mynowprog.eq.0.and.isb.eq.1)then                             11d29s17
         it=itrial(isk)                                                 12d23s17
         if(isk.eq.1)it=it+1                                            2d14s19
         if(idoub(isk)*iacto(isk).ne.0)then                             11d29s17
          do id=0,idoub(isk)-1                                           11d29s17
           iadb=iamat(isk)+iacto(isk)*id                                12d23s17
           do ia=0,iacto(isk)-1                                          11d29s17
            sumb=sumb+bc(iadb+ia)*bc(it+ia)                             12d23s17
           end do                                                       11d30s17
           it=it+iacto(isk)                                              11d29s17
          end do                                                         11d29s17
         end if                                                         11d29s17
        end if                                                          11d29s17
        kbmat=iamat(isk)+iacto(isk)*idoub(isk)                          12d23s17
        if(mynowprog.eq.0.and.isb.eq.1)then                             11d29s17
         if(noc(isk)*nvirtc(isk).ne.0)then                              11d29s17
          do io=0,noc(isk)-1                                            12d23s17
           iadb=kbmat+nvirtc(isk)*io                                    12d23s17
           do iv=0,nvirtc(isk)-1                                         11d29s17
            sumb=sumb+bc(iadb+iv)*bc(it+iv)                             12d23s17
           end do                                                       11d29s17
           it=it+nvirtc(isk)                                            12d23s17
          end do                                                        11d29s17
         end if                                                         11d29s17
        end if                                                          11d29s17
        iadmt=imt(isb)                                                  11d29s17
        if(isb.eq.1)iadmt=iadmt+1                                       11d29s17
c
c     hessa (a,b) idoub(b)idoub(a),iacto(a)iacto(b), a le b.
c
        it=itrial(isk)                                                  11d29s17
        if(isk.eq.1)it=it+1                                             11d29s17
        if(mynowprog.eq.0)then                                          11d29s17
         if(isk.le.isb)then                                             11d29s17
c
c     so a=isk and b=isb
c
          jhess=ihessa(isk,isb)                                         11d29s17
          do iab=0,iacto(isb)-1                                         11d29s17
           jadmt=iadmt+iab                                              11d29s17
           do iaa=0,iacto(isk)-1                                        11d29s17
            do ida=0,idoub(isk)-1                                       11d29s17
             jtrial=it+iaa+iacto(isk)*ida                               11d29s17
             fact=bc(jtrial)*xlami                                      11d29s17
             do idb=0,idoub(isb)-1                                      11d29s17
              iad=jadmt+iacto(isb)*idb
              bc(iad)=bc(iad)+bc(jhess+idb)*fact                        12d23s17
             end do                                                     11d29s17
             jhess=jhess+idoub(isb)                                     11d29s17
            end do                                                      11d29s17
           end do                                                       11d29s17
          end do                                                        11d29s17
         else                                                           11d29s17
c
c     so a=isb and b=isk
c
          do iab=0,iacto(isk)-1                                         11d29s17
           do iaa=0,iacto(isb)-1                                        11d29s17
            do ida=0,idoub(isb)-1                                       11d29s17
             do idb=0,idoub(isk)-1                                      11d29s17
              jtrial=it+iab+iacto(isk)*idb                              11d29s17
              kadmt=iadmt+iaa+iacto(isb)*ida                            11d29s17
              jhess=ihessa(isb,isk)+idb+idoub(isk)*(ida+idoub(isb)      11d30s17
     $             *(iaa+iacto(isb)*iab))                               11d29s17
              bc(kadmt)=bc(kadmt)+xlami*bc(jhess)*bc(jtrial)            11d29s17
             end do                                                     11d29s17
            end do                                                      11d29s17
           end do                                                       11d29s17
          end do                                                        11d29s17
         end if                                                         11d29s17
        end if                                                          11d29s17
c
c     hessb (a,b) act(a),doub(a);noc(b)nvirtc(b), all a and b.          11d30s17
c                                                                       11d29s17
        jadmt=iadmt+idoub(isb)*iacto(isb)                               11d29s17
        kadmt=imt(isk)                                                  11d29s17
        jt=itrial(isb)+iacto(isb)*idoub(isb)                            11d29s17
        if(isb.eq.1)jt=jt+1                                             11d29s17
        if(isk.eq.1)kadmt=kadmt+1                                       11d29s17
c
c     a=isk, b=isb
c
        call ilimts(noc(isb),nvirtc(isb),mynprocg,mynowprog,il,ih,      11d29s17
     $      i1s,i1e,i2s,i2e)                                            11d29s17
        i10=i1s                                                         11d29s17
        i1n=noc(isb)                                                    11d29s17
        nrow=iacto(isk)*idoub(isk)                                      11d29s17
        jhess=ihessb(isk,isb)                                           11d29s17
        do i2=i2s,i2e                                                   11d29s17
         if(i2.eq.i2e)i1n=i1e                                           11d29s17
         i2m=i2-1                                                       11d29s17
         do i1=i10,i1n                                                  11d29s17
          i1m=i1-1                                                      11d29s17
          id1=jadmt+i2m+nvirtc(isb)*i1m                                 11d29s17
          id2=jt+i2m+nvirtc(isb)*i1m                                    11d29s17
          fact=bc(id2)*xlami                                            11d29s17
          sum=0d0                                                       11d29s17
          do i34=0,nrow-1                                               11d29s17
           sum=sum+bc(it+i34)*bc(jhess+i34)                             11d29s17
           bc(kadmt+i34)=bc(kadmt+i34)+fact*bc(jhess+i34)               11d29s17
          end do                                                        11d29s17
          bc(id1)=bc(id1)+sum*xlami                                     11d29s17
          jhess=jhess+nrow                                              11d29s17
         end do                                                         11d29s17
         i10=1                                                          11d29s17
        end do                                                          11d29s17
c
c     hessc (a,b): noc(a),noc(b);nvirtc(b),nvirtc(isa), a le b          11d30s17
c     if a=b, we need to average in transpose as well
c
        it=it+iacto(isk)*idoub(isk)                                     11d30s17
        if(isk.lt.isb)then                                              11d30s17
c
c     a=isk, b=isb
c
         call ilimts(nvirtc(isb),nvirtc(isk),mynprocg,mynowprog,il,ih,  11d30s17
     $       i1s,i1e,i2s,i2e)                                           11d30s17
         i10=i1s                                                        11d30s17
         i1n=nvirtc(isb)                                                11d30s17
         nrow=noc(isb)*noc(isk)                                         11d30s17
         jhess=ihessc(isk,isb)                                          11d30s17
         do i2=i2s,i2e                                                  11d30s17
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e                                          11d30s17
          do i1=i10,i1n                                                 11d30s17
           i1m=i1-1
           do i4=0,noc(isb)-1                                           11d30s17
            id1=jadmt+i1m+nvirtc(isb)*i4                                11d30s17
            do i3=0,noc(isk)-1                                          11d30s17
             id2=it+i2m+nvirtc(isk)*i3                                  11d30s17
             id3=jhess+i3+noc(isk)*i4                                   11d30s17
             bc(id1)=bc(id1)+xlami*bc(id2)*bc(id3)                      12d28s17
  888        format(a1,4i5,f16.8,es15.7,2i5)
            end do                                                      11d30s17
           end do                                                       11d30s17
           jhess=jhess+nrow                                             11d30s17
          end do                                                        11d30s17
          i10=1                                                         11d30s17
         end do                                                         11d30s17
        else if(isb.lt.isk)then
c
c     a=isb, b=isk
c
         call ilimts(nvirtc(isk),nvirtc(isb),mynprocg,mynowprog,il,ih,  11d30s17
     $       i1s,i1e,i2s,i2e)                                           11d30s17
         i10=i1s                                                        11d30s17
         i1n=nvirtc(isk)                                                11d30s17
         nrow=noc(isb)*noc(isk)                                         11d30s17
         jhess=ihessc(isb,isk)                                          11d30s17
         do i2=i2s,i2e                                                  11d30s17
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e                                          11d30s17
          do i1=i10,i1n                                                 11d30s17
           i1m=i1-1
           do i4=0,noc(isk)-1                                           11d30s17
            do i3=0,noc(isb)-1                                          11d30s17
             id1=jadmt+i2m+nvirtc(isb)*i3                                11d30s17
             id2=it+i1m+nvirtc(isk)*i4                                  11d30s17
             id3=jhess+i3+noc(isb)*i4                                   11d30s17
             bc(id1)=bc(id1)+xlami*bc(id2)*bc(id3)                      12d28s17
            end do                                                      11d30s17
           end do                                                       11d30s17
           jhess=jhess+nrow                                             11d30s17
          end do                                                        11d30s17
          i10=1                                                         11d30s17
         end do                                                         11d30s17
        else
c
c     a=b=isb=isk
c
         call ilimts(nvirtc(isb),nvirtc(isb),mynprocg,mynowprog,il,ih,  11d30s17
     $       i1s,i1e,i2s,i2e)                                           11d30s17
         i10=i1s                                                        11d30s17
         i1n=nvirtc(isb)                                                11d30s17
         nrow=noc(isb)*noc(isb)                                         11d30s17
         jhess=ihessc(isb,isb)                                          11d30s17
         do i2=i2s,i2e                                                  11d30s17
          i2m=i2-1
          if(i2.eq.i2e)i1n=i1e                                          11d30s17
          do i1=i10,i1n                                                 11d30s17
           i1m=i1-1
           do i4=0,noc(isb)-1                                           11d30s17
            do i3=0,noc(isb)-1                                          11d30s17
             id1=jadmt+i1m+nvirtc(isb)*i4                                11d30s17
             id2=it+i2m+nvirtc(isb)*i3                                  11d30s17
             id3=jhess+i3+noc(isb)*i4                                   11d30s17
             jd1=jadmt+i2m+nvirtc(isb)*i3                                11d30s17
             jd2=it+i1m+nvirtc(isb)*i4                                  11d30s17
             bc(id1)=bc(id1)+0.5d0*xlami*bc(id2)*bc(id3)                11d30s17
             bc(jd1)=bc(jd1)+0.5d0*xlami*bc(jd2)*bc(id3)                11d30s17
            end do                                                      11d30s17
           end do                                                       11d30s17
           jhess=jhess+nrow                                             11d30s17
          end do                                                        11d30s17
          i10=1                                                         11d30s17
         end do                                                         11d30s17
        end if
       end do                                                           11d29s17
c
c     constraint times gradient
c
       if(mynowprog.eq.0)then                                           11d29s17
        if(isb.eq.1)bc(imt(1))=sumb                                     12d23s17
        iad=imt(isb)                                                    12d23s17
        if(isb.eq.1)iad=iad+1                                           2d14s19
        if(idoub(isb)*iacto(isb).ne.0)then                              11d29s17
         do id=0,idoub(isb)-1                                           11d29s17
          iadb=iamat(isb)+iacto(isb)*id                                 12d23s17
          do ia=0,iacto(isb)-1                                          11d29s17
           bc(iad+ia)=bc(iad+ia)+const*bc(iadb+ia)                      12d23s17
          end do                                                        11d29s17
          iad=iad+iacto(isb)                                            11d29s17
         end do                                                         11d29s17
        end if                                                          11d29s17
       end if                                                           11d29s17
       jbmat=iamat(isb)+iacto(isb)*idoub(isb)                           12d23s17
       if(mynowprog.eq.0)then                                           11d29s17
        if(noc(isb)*nvirtc(isb).ne.0)then                               11d29s17
         do io=0,noc(isb)-1                                             11d30s17
          do iv=0,nvirtc(isb)-1                                          11d29s17
           iadb=jbmat+iv+nvirtc(isb)*io                                 11d29s17
           bc(iad+iv)=bc(iad+iv)+const*bc(iadb)                         12d23s17
          end do                                                        11d29s17
          iad=iad+nvirtc(isb)                                           11d29s17
         end do                                                         11d29s17
        end if                                                          11d29s17
       end if                                                           11d29s17
      end do
      return
      end
