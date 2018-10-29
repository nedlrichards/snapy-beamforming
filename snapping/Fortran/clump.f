      subroutine clump(kbeam, vbeam, llasc, iclmp, vclmp, nbeam, nlabl)
      implicit none

      integer*4, intent(in)  :: kbeam(4, nbeam)
      integer*4, intent(in)  :: llasc(2, nbeam)
      integer*4, intent(in)  :: nbeam, nlabl
      real*8, intent(in)     :: vbeam(nbeam)

      integer*4, intent(out) :: iclmp(3, nlabl)
      real*8, intent(out)    :: vclmp(nlabl)

      integer*4 :: i, ivtmp, lcurr, ltrgt

      do 1 i=1, nlabl
          vclmp(i) = 0.
          iclmp(1, i) = 0
          iclmp(2, i) = 0
          iclmp(3, i) = 0
    1 continue

      do 2 i=1, nbeam
          lcurr = i
          ltrgt = llasc(1, lcurr)
C Find parent of current index
  801     if ( lcurr .NE. ltrgt ) then
              lcurr = ltrgt
              ltrgt = llasc(1, lcurr)
              go to 801
          end if
C Compare max value to current max
          ivtmp = llasc(2, lcurr)
          if ( vclmp(ivtmp) .LT. vbeam(i) ) then
              vclmp(ivtmp) = vbeam(i)
              iclmp(1, ivtmp) = kbeam(3, i)
              iclmp(2, ivtmp) = kbeam(4, i)
          end if
C account for number of elements in clump
          iclmp(3, ivtmp) = iclmp(3, ivtmp) + kbeam(2, i) - kbeam(1, i)
     ;                      + 1
    2 continue
      end
