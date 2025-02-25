      subroutine mkmem(maxbc)
      integer maxbc
      double precision, allocatable :: xmem(:)
      allocate (xmem(maxbc))
      call runit(xmem,xmem)
      call dws_synca
      deallocate (xmem)
      call dws_finalize
      stop
      end
