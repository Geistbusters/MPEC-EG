c mpec2.1 version eta. For copyright and Disclaimers, see start.f
      subroutine hcdd4v(vdoub,hcdoub,ncd,multh,iorb,nbasispc,natom,     11d27s18
     $     ngaus,ibdat,ipair,isym,iapair,ibstor,isstor,idorel,          4d24s20
     $     ascale,ndoub,nrootu,nbasisp,ntot,sr2,srh,bc,ibc)             11d9s22
      implicit real*8 (a-h,o-z)
c
c     superentend calculation of 4 virtual integral contribution to
c     hdd*cd
c
      include "common.store"
      include "common.mrci"
      include "common.hf"
      dimension multh(8,8),iorb(8),ncd(3,8,2),nbasispc(8),iaddr(36,2,4),  11d23s18
     $     nfcn(36,3),vdoub(*),hcdoub(*),iptoh(8,8),nbasisp(8)          11d10s20
      common/mycom/mynowprog,mynprocg,mynowprol,mynprol,mynownode,      5d24s18
     $     mynnode                                                      5d24s18
      ibcoffo=ibcoff
      call vdmo2so(vdoub,iorb,multh,ncd,nbasispc,iaddr,nfcn,nrootu,     4d9s19
     $     ntot,nbasisp,idorel,srh,0,bc,ibc)                            6d7s24
      iter=1
      idwsdeb=0
      isym=isymmrci                                                     12d7s18
      call paraerik(natom,ngaus,ibdat,idum,ihmat,ipair,nhcolt,isym,     4d24s20
     $     iapair,ibstor,isstor,multh,iptoh,iter,idwsdeb,idorel,ascale, 11d27s18
     $     nbasisp,iaddr,nfcn,nbasispc,ncd,iorb,hcdoub,ntot,nrootu,sr2, 6d22s21
     $     srh,bc,ibc)                                                  11d9s22
      return
      end
