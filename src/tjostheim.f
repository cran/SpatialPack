      subroutine tjostheim(xposx, xposy, yposx, yposy, n,
     &  xbar, ybar, cor)
c
c     arguments
      double precision  xposx(*), xposy(*), yposx(*), yposy(*)
      integer           n
      double precision  xbar, ybar, cor
c
c     local scalars
      double precision  sxx, syy, spx, spy
      integer           i
c
c     executable statements
c
      sxx = 0.d0
      syy = 0.d0
      spx = 0.d0
      spy = 0.d0
c
      do i = 1,n
        sxx = sxx + (xposx(i) - xbar) * (xposy(i) - xbar)
        syy = syy + (yposx(i) - ybar) * (yposy(i) - ybar)
        spx = spx + ((xposx(i) - xbar)**2 + (yposx(i) - ybar)**2)
        spy = spy + ((xposy(i) - xbar)**2 + (yposy(i) - ybar)**2)
      end do
c
      cor = (sxx + syy) / sqrt(spx * spy)
c
c     end of tjostheim
c
      end
