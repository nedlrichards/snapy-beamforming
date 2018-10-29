      subroutine rnext(kbeam, icurr, irnex, irmin, nbeam)
      implicit none
C Find the index of the first beam in the row next to the current one
C
C Used in the associate algorithms
C
C ::INPUTS::
C
C nbeams: total number of beams
C kbeams: index matrix of beam results
C icurr: beam number to start search on
C
      integer*4, intent(in)::nbeam, irmin
      integer*4, intent(in)::kbeam(4, nbeam), icurr
C
C ::OUTPUTS::
C
C irnex: index of first beam in the next row, 0 if we have exhausted beams.
C
      integer*4, intent(out)::irnex
C
      integer*4 :: i, ircur
Cf2py intent(hide), depend(nbeam) kbeam
      irnex = 0
C short circut if we are at the end of the beam output
      if ( icurr .GE. nbeam ) goto 299
C find the next time when row numbers dont match
      ircur = kbeam(4, icurr)
C loop untill the index changes or there are no more beams
      do 1 i = icurr + 1, nbeam
          if( kbeam(4, i) .NE. ircur .AND.
     ;        kbeam(4, i) .GE. irmin) then
              irnex = i
              go to 299
          end if
    1 continue
  299 end
