c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine xtimesn2(ncsfb,ncsfk,ncsfmid,mcol,iwpb,iwpk,xin,ndim,   12d4s20
     $     xout,ndim2,f1,f2,icsf0,mcsf,bc,ibc)                          11d10s22
c                                                                       2d22s21
c     if icsf0=1 and mcsf=ncsfk, this is the same as xtimesn.           2d22s21
c                                                                       2d22s21
      implicit real*8 (a-h,o-z)                                         7d17s20
      dimension xin(ndim,*),xout(ndim2,*)                               12d4s20
      include "common.store"                                            7d17s20
      if(mcol.le.0)return                                               11d3s22
      if(abs(f1).lt.1d-14)then                                          11d3s22
       if(f2.eq.0d0)then                                                11d3s22
        do i=1,mcol                                                      11d3s22
         do j=1,ncsfb                                                    11d3s22
          xout(j,i)=0d0                                                 11d3s22
         end do                                                           11d3s22
        end do                                                           11d3s22
       else                                                             11d3s22
        do i=1,mcol                                                      11d3s22
         do j=1,ncsfb                                                    11d3s22
          xout(j,i)=xout(j,i)*f2                                        11d3s22
         end do                                                           11d3s22
        end do                                                           11d3s22
       end if                                                           11d3s22
       return                                                           11d3s22
      end if                                                            11d3s22
      ibcoffo=ibcoff                                                    11d3s22
      n1=ncsfmid*ncsfk*mcol                                             11d3s22
      if(iwpk.lt.0)n1=n1/2                                              11d3s22
      n2=ncsfb*ncsfmid*mcol                                             11d3s22
      if(iwpb.lt.0)n2=n2/2                                              11d3s22
      nop1=n1+n2                                                        11d3s22
      n1=ncsfb*ncsfmid*ncsfk                                            11d3s22
      if(iwpb.lt.0)n1=n1/2                                              11d3s22
      if(iwpk.lt.0)n1=n1*3                                              11d3s22
      nop2=n1+ncsfb*ncsfk*mcol                                          11d3s22
      if(nop1.lt.nop2)then                                              11d3s22
       itmp1=ibcoff                                                     12d5s20
       ibcoff=itmp1+ncsfmid*mcol                                        12d5s20
       if(iwpk.lt.0)then                                                12d5s20
        nusedi=ibc(-iwpk)/2                                             12d4s20
        if(2*nusedi.ne.ibc(-iwpk))nusedi=nusedi+1                       12d4s20
        icmp1=-iwpk+1                                                   12d4s20
        icmp2=icmp1+ncsfk                                               12d5s20
        icmp3=icmp2+nusedi                                              12d4s20
        call compxtimes2(ibc(icmp1),ibc(icmp2),bc(icmp3),xin,            12d5s20
     $       ndim,ncsfk,bc(itmp1),ncsfmid,ncsfmid,mcol,0d0,icsf0,mcsf,  11d14s22
     $       bc,ibc)                                                    11d14s22
        do i=0,ncsfmid*mcol-1                                           12d5s20
         bc(itmp1+i)=bc(itmp1+i)*f1                                     12d5s20
        end do                                                          12d5s20
       else                                                             12d5s20
        jwpk=iwpk+ncsfmid*(icsf0-1)                                     2d22s21
        call dgemm('n','n',ncsfmid,mcol,mcsf,f1,                        2d22s21
     $            bc(jwpk),ncsfmid,xin,ndim,0d0,                        12d5s20
     $            bc(itmp1),ncsfmid,                                    12d5s20
     d' xtimesn2.  1')
       end if                                                           12d5s20
       if(iwpb.lt.0)then                                                12d5s20
        nusedi=ibc(-iwpb)/2                                             12d4s20
        if(2*nusedi.ne.ibc(-iwpb))nusedi=nusedi+1                       12d4s20
        icmp1=-iwpb+1                                                   12d4s20
        icmp2=icmp1+ncsfmid                                             12d5s20
        icmp3=icmp2+nusedi                                              12d4s20
        call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itmp1),      12d5s20
     $       ncsfmid,ncsfmid,xout,ndim2,ncsfb,mcol,f2,1d0)              11d3s22
       else                                                             12d5s20
        call dgemm('n','n',ncsfb,mcol,ncsfmid,1d0,                      12d5s20
     $      bc(iwpb),ncsfb,bc(itmp1),ncsfmid,f2,xout,ndim2,             11d3s22
     d' xtimesn2.  2')
       end if                                                           12d5s20
       ibcoff=itmp1                                                     12d5s20
      else                                                              11d3s22
       itmp1=ibcoff                                                     11d3s22
       ibcoff=itmp1+ncsfb*mcsf                                          11d3s22
       ibcoff=itmp1+ncsfb*ncsfk                                         11d3s22
       call enough('xtimes2.3',bc,ibc)
       if(iwpk.lt.0)then                                                11d3s22
        itrial=ibcoff                                                   11d3s22
        ibcoff=itrial+ncsfmid*ncsfk                                     11d3s22
        call enough('xtimes2.4',bc,ibc)
        nusedi=ibc(-iwpk)/2                                             12d4s20
        if(2*nusedi.ne.ibc(-iwpk))nusedi=nusedi+1                       12d4s20
        icmp1=-iwpk+1                                                   12d4s20
        icmp2=icmp1+ncsfk                                               12d5s20
        icmp3=icmp2+nusedi                                              12d4s20
        call uncompxu(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(itrial),       9d30s22
     $        ncsfmid,ncsfk)                                            11d3s22
       else                                                             11d3s22
        itrial=iwpk                                                     11d3s22
       end if                                                           11d3s22
       jtrial=itrial+ncsfmid*(icsf0-1)                                  11d3s22
       if(iwpb.lt.0)then                                                12d5s20
        nusedi=ibc(-iwpb)/2                                             12d4s20
        if(2*nusedi.ne.ibc(-iwpb))nusedi=nusedi+1                       12d4s20
        icmp1=-iwpb+1                                                   12d4s20
        icmp2=icmp1+ncsfmid                                             12d5s20
        icmp3=icmp2+nusedi                                              12d4s20
        call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(jtrial),      12d5s20
     $       ncsfmid,ncsfmid,bc(itmp1),ncsfb,ncsfb,mcsf,0d0,1d0)              3d25s22
       else                                                             12d5s20
        call dgemm('n','n',ncsfb,mcsf,ncsfmid,1d0,                      11d3s22
     $      bc(iwpb),ncsfb,bc(jtrial),ncsfmid,0d0,bc(itmp1),ncsfb,      11d3s22
     d' xtimesn2.  3')                                                  11d3s22
       end if                                                           12d5s20
       call dgemm('n','n',ncsfb,mcol,mcsf,f1,bc(itmp1),ncsfb,           11d3s22
     $      xin,ndim,f2,xout,ndim2,                                     11d3s22
     $      'xtimesn2. 4')                                              11d3s22
       if(iwpk.lt.0)ibcoff=itrial                                       11d3s22
      end if                                                            11d3s22
      ibcoff=ibcoffo                                                    11d3s22
      return                                                            7d17s20
      end                                                               7d17s20
