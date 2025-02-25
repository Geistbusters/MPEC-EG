c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine mtimesx(z,p,ihessa,ihessb,ihessc,nsymb,multh,ipuse,    5d9s22
     $     ipt,idoub,iacto,noc,nvirt,nvar,bc,ibc)                       11d28s22
      implicit real*8 (a-h,o-z)                                         5d9s22
      dimension z(*),p(*),ihessa(8,*),ihessb(8,*),ihessc(8,*),          5d9s22
     $     multh(8,8),ipt(8,2),idoub(*),iacto(*),noc(*),nvirt(*)        11d28s22
      include "common.store"                                            5d9s22
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,
     $     mynnode
      do i=1,nvar                                                       5d9s22
       z(i)=0d0
      end do
      do isk=1,nsymb                                                    5d9s22
       isko=multh(isk,ipuse)                                            5d9s22
       do isb=1,nsymb                                                   5d9s22
        isbo=multh(isb,ipuse)                                           5d9s22
        if(ihessa(isb,isk).gt.0.and.mynowprog.eq.0)then                 5d9s22
c     hessa (a,b) idoub(b)idoub(a),iacto(a)iacto(b), a le b.
         ii=ihessa(isb,isk)                                             5d9s22
         jp=1+ipt(isk,1)                                                5d9s22
         jz=1+ipt(isb,1)                                                5d9s22
         do i4=0,iacto(isk)-1                                           5d9s22
          do i3=0,iacto(isb)-1                                          5d9s22
           do i2=0,idoub(isbo)-1                                        5d9s22
            do i1=0,idoub(isko)-1                                       5d9s22
             iadp=jp+i4+iacto(isk)*i1                                   5d9s22
             iadz=jz+i3+iacto(isb)*i2                                   5d9s22
             z(iadz)=z(iadz)+bc(ii)*p(iadp)                             5d9s22
             ii=ii+1                                                    5d9s22
            end do                                                      5d9s22
           end do                                                       5d9s22
          end do                                                        5d9s22
         end do                                                         5d9s22
         if(isb.ne.isk)then                                             5d9s22
          ii=ihessa(isb,isk)                                            5d9s22
          jp=1+ipt(isb,1)                                               5d9s22
          jz=1+ipt(isk,1)                                               5d9s22
          do i4=0,iacto(isk)-1                                          5d9s22
           do i3=0,iacto(isb)-1                                         5d9s22
            do i2=0,idoub(isbo)-1                                       5d9s22
             do i1=0,idoub(isko)-1                                      5d9s22
              iadp=jp+i3+iacto(isb)*i2                                  5d9s22
              iadz=jz+i4+iacto(isk)*i1                                  5d9s22
              z(iadz)=z(iadz)+bc(ii)*p(iadp)                            5d9s22
              ii=ii+1                                                   5d9s22
             end do                                                     5d9s22
            end do                                                      5d9s22
           end do                                                       5d9s22
          end do                                                        5d9s22
         end if                                                         5d9s22
        end if                                                          5d9s22
        if(ihessb(isb,isk).gt.0)then                                    5d9s22
c     hessb (a,b) act(a),doub(a);noc(b)nvirtc(b), all a and b.          11d30s17
         nrowb=idoub(isbo)*iacto(isb)                                   5d9s22
         call ilimts(noc(isko),nvirt(isk),mynprocg,mynowprog,ilb,ihb,   5d9s22
     $         i1sb,i1eb,i2sb,i2eb)
         nhereb=ihb+1-ilb                                               5d9s22
         if(min(nrowb,nhereb).gt.0)then                                 5d9s22
          i10=i1sb                                                      5d9s22
          i1n=noc(isko)                                                 5d9s22
          ii=ihessb(isb,isk)                                            5d9s22
          jp=1+ipt(isk,2)                                               5d9s22
          jz=1+ipt(isb,1)                                               5d9s22
          do i2=i2sb,i2eb                                               5d9s22
           i2m=i2-1                                                     5d9s22
           if(i2.eq.i2eb)i1n=i1eb                                       5d9s22
           do i1=i10,i1n                                                5d9s22
            i1m=i1-1                                                    5d9s22
            do i4=0,idoub(isbo)-1                                       5d9s22
             do i3=0,iacto(isb)-1                                       5d9s22
              iadp=jp+i2m+nvirt(isk)*i1m                                5d13s22
              iadz=jz+i3+iacto(isb)*i4                                  5d9s22
              z(iadz)=z(iadz)+bc(ii)*p(iadp)                            5d9s22
              ii=ii+1                                                   5d9s22
             end do                                                     5d9s22
            end do                                                      5d9s22
           end do                                                       5d9s22
           i10=1                                                        5d9s22
          end do                                                        5d9s22
          i10=i1sb                                                      5d9s22
          i1n=noc(isko)                                                 5d9s22
          ii=ihessb(isb,isk)                                            5d9s22
          jp=1+ipt(isb,1)                                               5d9s22
          jz=1+ipt(isk,2)                                               5d9s22
          do i2=i2sb,i2eb                                               5d9s22
           i2m=i2-1                                                     5d9s22
           if(i2.eq.i2eb)i1n=i1eb                                       5d9s22
           do i1=i10,i1n                                                5d9s22
            i1m=i1-1                                                    5d9s22
            do i4=0,idoub(isbo)-1                                       5d9s22
             do i3=0,iacto(isb)-1                                       5d9s22
              iadp=jp+i3+iacto(isb)*i4                                  5d9s22
              iadz=jz+i2m+nvirt(isk)*i1m                                5d13s22
              z(iadz)=z(iadz)+bc(ii)*p(iadp)                            5d9s22
              ii=ii+1                                                   5d9s22
             end do                                                     5d9s22
            end do                                                      5d9s22
           end do                                                       5d9s22
           i10=1                                                        5d9s22
          end do                                                        5d9s22
         end if                                                         5d9s22
        end if                                                          5d9s22
        if(ihessc(isb,isk).gt.0)then                                    5d9s22
c
c     to avoid extra communication, hessc is not symmetric, however
c     0.5*(hessc+hesscT) is the correct result.
c     hessc (a,b): noc(a),noc(b);nvirtc(b),nvirtc(isa), a le b          11d30s17
c     if a=b, we need to average in transpose as well
         nrowc=noc(isbo)*noc(isko)                                      5d9s22
         call ilimts(nvirt(isk),nvirt(isb),mynprocg,mynowprog,ilh,ihh,  5d9s22
     $         i1sh,i1eh,i2sh,i2eh)                                     5d9s22
         nherec=ihh+1-ilh                                               5d9s22
         if(min(nrowc,nherec).gt.0)then                                 5d9s22
          ii=ihessc(isb,isk)                                            5d9s22
          i10=i1sh                                                      5d9s22
          i1n=nvirt(isk)                                                5d9s22
          jp=1+ipt(isk,2)                                               5d9s22
          jz=1+ipt(isb,2)                                               5d9s22
          if(isb.eq.isk)then                                            5d13s22
           fmul=0.5d0                                                   5d13s22
          else                                                          5d13s22
           fmul=1d0                                                     5d13s22
          end if                                                        5d13s22
          do i2=i2sh,i2eh                                               5d9s22
           i2m=i2-1                                                     5d9s22
           if(i2.eq.i2eh)i1n=i1eh                                       5d9s22
           do i1=i10,i1n                                                5d9s22
            i1m=i1-1                                                    5d9s22
            do i4=0,noc(isko)-1                                         5d9s22
             do i3=0,noc(isbo)-1                                        5d9s22
              iadp=jp+i1m+nvirt(isk)*i4                                 5d9s22
              iadz=jz+i2m+nvirt(isb)*i3                                 5d9s22
              z(iadz)=z(iadz)+fmul*bc(ii)*p(iadp)                       5d13s22
              ii=ii+1                                                   5d9s22
             end do                                                     5d9s22
            end do                                                      5d9s22
           end do                                                       5d9s22
           i10=1                                                        5d9s22
          end do                                                        5d9s22
           ii=ihessc(isb,isk)                                           5d9s22
           i10=i1sh                                                     5d9s22
           i1n=nvirt(isk)                                               5d9s22
           jp=1+ipt(isb,2)                                              5d9s22
           jz=1+ipt(isk,2)                                              5d9s22
           do i2=i2sh,i2eh                                              5d9s22
            i2m=i2-1                                                    5d9s22
            if(i2.eq.i2eh)i1n=i1eh                                      5d9s22
            do i1=i10,i1n                                               5d9s22
             i1m=i1-1                                                   5d9s22
             do i4=0,noc(isko)-1                                        5d9s22
              do i3=0,noc(isbo)-1                                       5d9s22
               iadp=jp+i2m+nvirt(isb)*i3                                5d9s22
               iadz=jz+i1m+nvirt(isk)*i4                                5d9s22
               z(iadz)=z(iadz)+fmul*bc(ii)*p(iadp)                           5d9s22
               ii=ii+1                                                  5d9s22
              end do                                                    5d9s22
             end do                                                     5d9s22
            end do                                                      5d9s22
            i10=1                                                       5d9s22
           end do                                                       5d9s22
         end if                                                         5d9s22
        end if                                                          5d9s22
       end do                                                           5d9s22
      end do                                                            5d9s22
      call dws_gsumf(z,nvar)
      return                                                            5d9s22
      end                                                               5d9s22
