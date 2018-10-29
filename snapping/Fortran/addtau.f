      subroutine addtau(ichan, itau, ival, mshif, nchan, nnz)
      implicit none
C Apply a constant shift to index vector by channel
      integer*4, intent(in)::nnz, nchan
      integer*4, intent(in)::ichan(nchan + 1), itau(nchan), ival(nnz)
      integer*4, intent(out)::mshif(nnz)
      integer*4 :: i, j
Cf2py intent(hide), depend(nchan + 1) ichan
Cf2py intent(hide), depend(nchan) itau
Cf2py intent(hide), depend(nnz) ival
Cf2py intent(hide), depend(nnz) mshif

C Add channel delays to each channel and save to mshif
      do 1 i = 1, nchan
          do 1 j = ichan(i), (ichan(i + 1) - 1)
              mshif(j) = ival(j) - itau(i)
    1 continue
      end
