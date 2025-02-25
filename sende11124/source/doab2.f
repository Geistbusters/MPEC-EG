c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine doab2(idatac,mca,idatbc,mcb,nadetb,nadetk,nbdetb,      3d16s17
     $     nbdetk,nroot,gb,nhereb,ilb,ihb,gk,nherek,ilk,ihk,vecab,vecak,3d16s17
     $     vecbb,vecbk,i2e,irelo,ism,iacto,iu1,iu2,mnza,mnzb,nza,nzb,   11d14s22
     $     bc,ibc)                                                      11d14s22
      implicit real*8 (a-h,o-z)                                         3d16s17
      integer*2 idatac(4,mca),idatbc(4,mcb),idwstmp2                    3d16s17
      integer*1 idwstmp1(4)                                             3d16s17
      equivalence (idwstmp2,idwstmp1)                                   3d16s17
      integer*8 ism(*),irelo(*)                                         3d16s17
      include "common.store"                                            3d16s17
      logical log1,log2                                                 4d9s18
      dimension gb(nroot,nhereb,nbdetb),gk(nroot,nherek,nbdetk),i2e(*), 3d16s17
     $     iacto(8),phs(2),vecab(nroot,nadetb,nbdetb),                  3d16s17
     $     vecak(nroot,nadetk,nbdetk),mnza(3,*),mnzb(3,*)               4d30s18
      data phs/1d0,-1d0/                                                3d16s17
      common/fnd2cm/inv(2,8,8,8)                                        4d9s18
      ioa=0                                                             4d30s18
      igoal=2058958
      do ita=1,nza                                                      4d30s18
       iab=mnza(1,ita)                                                  4d30s18
       nrow=(iacto(iab)*(iacto(iab)+1))/2                               4d30s18
       iak=mnza(2,ita)                                                  4d30s18
       iob=0
       do itb=1,nzb                                                     4d30s18
        ibb=mnzb(iu1,itb)
        ibk=mnzb(iu2,itb)                                               4d30s18
        i2eu=inv(1,iab,iak,ibb)                                         4d30s18
        icase=inv(2,iab,iak,ibb)                                        4d30s18
        if(iab.eq.iak)then                                              4d30s18
         do ia0=1,mnza(3,ita)                                             4d30s18
          ia=ia0+ioa                                                      4d30s18
          i1a=idatac(1,ia)
          i2a=idatac(2,ia)
          log1=i1a.ge.ilb.and.i1a.le.ihb
          fact1=0d0                                                     4d30s18
          if(log1)fact1=1d0                                             4d30s18
          log2=i2a.ge.ilk.and.i2a.le.ihk
          fact2=0d0                                                     4d30s18
          if(log2)fact2=1d0                                             4d30s18
          if(log1.or.log2)then                                             4d9s18
           idwstmp2=idatac(3,ia)
           nab=irelo(idwstmp1(1))
           nak=irelo(idwstmp1(2))
           pa=phs(idatac(4,ia))
           i2ao=i2a-ilk+1
           i1ao=i1a-ilb+1
           ix=max(nab,nak)
           in=min(nab,nak)
           ii1=((ix*(ix-1))/2)+in-1                                      3d16s17
           iad1=i2e(i2eu)+ii1
           if(i1ao.ge.1.and.i1ao.le.nhereb)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             ix=max(nbb,nbk)
             in=min(nbb,nbk)
             ii2=iad1+nrow*(((ix*(ix-1))/2)+in-1)                        4d30s18
             fact=pb*bc(ii2)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gb(k,i1ao,i1b)=gb(k,i1ao,i1b)+facta*vecak(k,i2a,i2b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
           if(i2ao.ge.1.and.i2ao.le.nherek)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             ix=max(nbb,nbk)
             in=min(nbb,nbk)
             ii2=iad1+nrow*(((ix*(ix-1))/2)+in-1)                        4d30s18
             fact=pb*bc(ii2)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gk(k,i2ao,i2b)=gk(k,i2ao,i2b)+factb*vecab(k,i1a,i1b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
          end if                                                        4d30s18
         end do                                                         4d30s18
        else if(icase.eq.1)then                                         4d30s18
         do ia0=1,mnza(3,ita)                                             4d30s18
          ia=ia0+ioa                                                      4d30s18
          i1a=idatac(1,ia)
          i2a=idatac(2,ia)
          log1=i1a.ge.ilb.and.i1a.le.ihb
          fact1=0d0                                                     4d30s18
          if(log1)fact1=1d0                                             4d30s18
          log2=i2a.ge.ilk.and.i2a.le.ihk
          fact2=0d0                                                     4d30s18
          if(log2)fact2=1d0                                             4d30s18
          if(log1.or.log2)then                                             4d9s18
           idwstmp2=idatac(3,ia)
           nab=irelo(idwstmp1(1))
           nak=irelo(idwstmp1(2))
           pa=phs(idatac(4,ia))
           i2ao=i2a-ilk+1
           i1ao=i1a-ilb+1
           if(i1ao.ge.1.and.i1ao.le.nhereb)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             iad1=i2e(i2eu)+nab-1+iacto(iab)*(nak-1+iacto(iak)*(nbb-1     3d16s17
     $          +iacto(ibb)*(nbk-1)))                                    3d16s17
             fact=pb*bc(iad1)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gb(k,i1ao,i1b)=gb(k,i1ao,i1b)+facta*vecak(k,i2a,i2b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
           if(i2ao.ge.1.and.i2ao.le.nherek)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             iad1=i2e(i2eu)+nab-1+iacto(iab)*(nak-1+iacto(iak)*(nbb-1     3d16s17
     $          +iacto(ibb)*(nbk-1)))                                    3d16s17
             fact=pb*bc(iad1)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gk(k,i2ao,i2b)=gk(k,i2ao,i2b)+factb*vecab(k,i1a,i1b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
          end if                                                        4d30s18
         end do                                                         4d30s18
        else if(icase.eq.2)then                                         4d30s18
         do ia0=1,mnza(3,ita)                                             4d30s18
          ia=ia0+ioa                                                      4d30s18
          i1a=idatac(1,ia)
          i2a=idatac(2,ia)
          log1=i1a.ge.ilb.and.i1a.le.ihb
          fact1=0d0                                                     4d30s18
          if(log1)fact1=1d0                                             4d30s18
          log2=i2a.ge.ilk.and.i2a.le.ihk
          fact2=0d0                                                     4d30s18
          if(log2)fact2=1d0                                             4d30s18
          if(log1.or.log2)then                                             4d9s18
           idwstmp2=idatac(3,ia)
           nab=irelo(idwstmp1(1))
           nak=irelo(idwstmp1(2))
           pa=phs(idatac(4,ia))
           i2ao=i2a-ilk+1
           i1ao=i1a-ilb+1
           if(i1ao.ge.1.and.i1ao.le.nhereb)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             iad1=i2e(i2eu)+nak-1+iacto(iak)*(nab-1+iacto(iab)*(nbk-1     3d16s17
     $          +iacto(ibk)*(nbb-1)))                                    3d16s17
             fact=pb*bc(iad1)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gb(k,i1ao,i1b)=gb(k,i1ao,i1b)+facta*vecak(k,i2a,i2b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
           if(i2ao.ge.1.and.i2ao.le.nherek)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             iad1=i2e(i2eu)+nak-1+iacto(iak)*(nab-1+iacto(iab)*(nbk-1     3d16s17
     $          +iacto(ibk)*(nbb-1)))                                    3d16s17
             fact=pb*bc(iad1)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gk(k,i2ao,i2b)=gk(k,i2ao,i2b)+factb*vecab(k,i1a,i1b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
          end if                                                        4d30s18
         end do                                                         4d30s18
        else if(icase.eq.3)then                                         4d30s18
         do ia0=1,mnza(3,ita)                                             4d30s18
          ia=ia0+ioa                                                      4d30s18
          i1a=idatac(1,ia)
          i2a=idatac(2,ia)
          log1=i1a.ge.ilb.and.i1a.le.ihb
          fact1=0d0                                                     4d30s18
          if(log1)fact1=1d0                                             4d30s18
          log2=i2a.ge.ilk.and.i2a.le.ihk
          fact2=0d0                                                     4d30s18
          if(log2)fact2=1d0                                             4d30s18
          if(log1.or.log2)then                                             4d9s18
           idwstmp2=idatac(3,ia)
           nab=irelo(idwstmp1(1))
           nak=irelo(idwstmp1(2))
           pa=phs(idatac(4,ia))
           i2ao=i2a-ilk+1
           i1ao=i1a-ilb+1
           if(i1ao.ge.1.and.i1ao.le.nhereb)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             iad1=i2e(i2eu)+nab-1+iacto(iab)*(nak-1+iacto(iak)*(nbk-1    4d30s18
     $          +iacto(ibk)*(nbb-1)))                                    3d16s17
             fact=pb*bc(iad1)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gb(k,i1ao,i1b)=gb(k,i1ao,i1b)+facta*vecak(k,i2a,i2b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
           if(i2ao.ge.1.and.i2ao.le.nherek)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             iad1=i2e(i2eu)+nab-1+iacto(iab)*(nak-1+iacto(iak)*(nbk-1    4d30s18
     $          +iacto(ibk)*(nbb-1)))                                    3d16s17
             fact=pb*bc(iad1)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gk(k,i2ao,i2b)=gk(k,i2ao,i2b)+factb*vecab(k,i1a,i1b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
          end if                                                        4d30s18
         end do                                                         4d30s18
        else                                                            4d30s18
         do ia0=1,mnza(3,ita)                                             4d30s18
          ia=ia0+ioa                                                      4d30s18
          i1a=idatac(1,ia)
          i2a=idatac(2,ia)
          log1=i1a.ge.ilb.and.i1a.le.ihb
          fact1=0d0                                                     4d30s18
          if(log1)fact1=1d0                                             4d30s18
          log2=i2a.ge.ilk.and.i2a.le.ihk
          fact2=0d0                                                     4d30s18
          if(log2)fact2=1d0                                             4d30s18
          if(log1.or.log2)then                                             4d9s18
           idwstmp2=idatac(3,ia)
           nab=irelo(idwstmp1(1))
           nak=irelo(idwstmp1(2))
           pa=phs(idatac(4,ia))
           i2ao=i2a-ilk+1
           i1ao=i1a-ilb+1
           if(i1ao.ge.1.and.i1ao.le.nhereb)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             iad1=i2e(i2eu)+nak-1+iacto(iak)*(nab-1+iacto(iab)*(nbb-1     3d16s17
     $          +iacto(ibb)*(nbk-1)))                                    3d16s17
             fact=pb*bc(iad1)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gb(k,i1ao,i1b)=gb(k,i1ao,i1b)+facta*vecak(k,i2a,i2b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
           if(i2ao.ge.1.and.i2ao.le.nherek)then
            do ib0=1,mnzb(3,itb)                                          4d30s18
             ib=ib0+iob                                                   4d30s18
             i1b=idatbc(iu1,ib)
             i2b=idatbc(iu2,ib)
             idwstmp2=idatbc(3,ib)
             nbb=irelo(idwstmp1(iu1))
             nbk=irelo(idwstmp1(iu2))
             pb=pa*phs(idatbc(4,ib))
             iad1=i2e(i2eu)+nak-1+iacto(iak)*(nab-1+iacto(iab)*(nbb-1     3d16s17
     $          +iacto(ibb)*(nbk-1)))                                    3d16s17
             fact=pb*bc(iad1)                                             4d30s18
             facta=fact1*fact                                            4d30s18
             factb=fact2*fact                                            4d30s18
             do k=1,nroot                                                4d30s18
              gk(k,i2ao,i2b)=gk(k,i2ao,i2b)+factb*vecab(k,i1a,i1b)         4d30s18
             end do                                                      4d30s18
            end do                                                        4d30s18
           end if
          end if                                                        4d30s18
         end do                                                         4d30s18
        end if                                                          4d30s18
        iob=iob+mnzb(3,itb)
       end do                                                           4d30s18
       ioa=ioa+mnza(3,ita)                                              4d30s18
      end do                                                            4d30s18
      return
      end
