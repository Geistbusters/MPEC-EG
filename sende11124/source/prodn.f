c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine prodn(iwpb,iwpk,ncsfb,ncsfk,ncsfmid,prod,bc,ibc,fi,fo) 2d13s23
      implicit real*8 (a-h,o-z)                                         12d4s20
      dimension prod(*)                                                 2d9s23
      include "common.store"                                            12d4s20
      ibctop=ibcoff                                                     12d5s20
      if(min(ncsfb,ncsfk,ncsfmid).le.0)then                             4d14s21
       write(6,*)('bad input to prodn: '),ncsfb,ncsfk,ncsfmid
       stop 'prodn'
      end if                                                            4d14s21
      iprod=ibcoff                                                      2d10s23
      itmp=iprod+ncsfb*ncsfk                                            2d10s23
      ibcoff=itmp+ncsfmid                                               2d10s23
      call enough('prodn.z',bc,ibc)                                     2d10s23
      if(iwpk.lt.0)then                                                 2d10s23
       nusedi=ibc(-iwpk)/2                                              12d5s20
       if(2*nusedi.ne.ibc(-iwpk))nusedi=nusedi+1                        12d5s20
       icmp1=-iwpk+1                                                    12d5s20
       icmp2=icmp1+ncsfk                                                12d5s20
       icmp3=icmp2+nusedi                                               12d4s20
       if(iwpb.lt.0)then                                                2d10s23
        nusedi=ibc(-iwpb)/2                                              12d5s20
        if(2*nusedi.ne.ibc(-iwpb))nusedi=nusedi+1                        12d5s20
        icmp1a=-iwpb+1                                                    12d5s20
        icmp2a=icmp1a+ncsfmid                                              12d5s20
        icmp3a=icmp2a+nusedi                                               12d5s20
        call compcomp(ncsfb,ibc(icmp1a),ibc(icmp2a),bc(icmp3a),ncsfmid,  2d9s23
     $      ibc(icmp1),ibc(icmp2),bc(icmp3),prod,ncsfb,ncsfk,fo,        2d13s23
     $      fi,bc(itmp))                                                2d13s23
       else                                                             2d10s23
        call timescomp(bc(iwpb),ncsfb,ncsfb,ibc(icmp1),ibc(icmp2),       2d9s23
     $      bc(icmp3),prod,ncsfb,ncsfk,fo,fi)                           2d13s23
       end if                                                           2d10s23
      else                                                              2d10s23
       if(iwpb.lt.0)then                                                2d10s23
        nusedi=ibc(-iwpb)/2                                              12d5s20
        if(2*nusedi.ne.ibc(-iwpb))nusedi=nusedi+1                        12d5s20
        icmp1=-iwpb+1                                                    12d5s20
        icmp2=icmp1+ncsfmid                                              12d5s20
        icmp3=icmp2+nusedi                                               12d5s20
        call compxtimes(ibc(icmp1),ibc(icmp2),bc(icmp3),bc(iwpk),       2d10s23
     $       ncsfmid,ncsfmid,prod,ncsfb,ncsfb,ncsfk,fo,fi)              2d13s23
       else                                                              12d5s20
        call dgemm('n','n',ncsfb,ncsfk,ncsfmid,fi,bc(iwpb),ncsfb,       2d13s23
     $     bc(iwpk),ncsfmid,fo,prod,ncsfb,                              2d13s23
     d' prodn.  1')
       end if                                                           2d10s23
      end if                                                            2d10s23
      ibcoff=ibctop                                                     12d5s20
      return                                                            12d4s20
      end                                                               12d4s20
