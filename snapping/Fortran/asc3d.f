      subroutine asc3d(kbeam, llasc, nrows, nbeam, nclmp)
      implicit none
C Associate beams which are next to each other
C
C Elements which are next to each other in any of three dimensions are
C associated into the same cluster. This is a connected component labeling
C algorithm with 8 point connectivity (generalized to 3D).
C
C ::INPUTS::
C
C kbeam: beamformer output of clump routine
C vval: nnz vector of beamformer values
C ichan: nchan vector marking the start of each channels index and value vector
C nnz: total number of index and value pairs
C nchan: number of data channels
C thresh: beamformer threshold. Only save values which sum GE this value
C lnmax: total allocated beamformer result positions
C

      integer*4, intent(in) :: nrows, nbeam
      integer*4, intent(in) :: kbeam(4, nbeam)
C
C ::OUTPUTS::
C llasc: 2 x nbeam matrix of label tree
C        (parent index, parent label number)
C nclmp: number of clusters
C
      integer*4, intent(out) :: llasc(2, nbeam), nclmp
C
C ::Internal::
C
      integer*4 :: i,krlfc, krupc, lcurr, kcurr, knext
      integer*4 :: ircur, irlfc, irnex, irupc, iaftr

C preallocate association array
      do 1 i=1, nbeam
          llasc(1, i) = i
          llasc(2, i) = 0
    1 continue

C start loop at first index
      ircur = 1
      kcurr = kbeam(4, 1)
      krupc = 0
      krlfc = 0
      irlfc = 0

C find the first index of the next row
      call rnext(kbeam, ircur, irnex, 0, nbeam)
C find the first index of the first row in next column
      call rnext(kbeam, ircur, irupc, nrows, nbeam)
      if ( irupc .GT. 0 ) then
          krupc = kbeam(4, irupc)
          irlfc = irupc
          krlfc = krupc
      end if

C 0 indicates we are at the end of the beams
  801 if ( irnex .LE. 0 ) go to 290
C make sure left up index is not too small
      if ( krlfc .LT. kcurr + nrows - 1 ) then
          call rnext(kbeam, irlfc, irlfc, kcurr + nrows - 1, nbeam)
          if ( irlfc .GT. 0 ) krlfc = kbeam(4, irlfc)
      end if
C step left up index by at least 1 if we are at the end of a row
  855 if ( irlfc .GT. 0 .AND. mod(krlfc, nrows) .EQ. nrows ) then
          call rnext(kbeam, irlfc, irlfc, 0, nbeam)
          if ( irlfc .GT. 0 ) krlfc = kbeam(4, irlfc)
          go to 855
      end if
C Make sure that up is not too small
      if ( krupc .LT. kcurr + nrows ) then
          call rnext(kbeam, ircur, irupc, kcurr + nrows, nbeam)
          if ( irupc .GT. 0 ) krupc = kbeam(4, irupc)
      end if
C associate left up if that is an option
      if ( krlfc .EQ. kcurr + nrows - 1 ) then
          call ascnx(kbeam, llasc, nbeam, nrows - 1,
     ;               ircur, irlfc, iaftr)
          irlfc = iaftr
          if ( irlfc .GT. 0 ) krlfc = kbeam(4, irlfc)
      end if
C see if one row up is an option
      if ( krupc .EQ. kcurr + nrows ) then
          call ascnx(kbeam, llasc, nbeam, nrows, ircur, irupc, iaftr)
          irupc = iaftr
          if ( irupc .GT. 0 ) krupc = kbeam(4, irupc)
      end if
C see if one row up is an option (right diagonal)
C need to check additionally if diagonal is in the current row
      if ( krupc .EQ. kcurr + nrows + 1
     ;     .AND. mod(krupc, nrows) .NE. 1 ) then
          call ascnx(kbeam, llasc, nbeam, nrows + 1,
     ;               ircur, irupc, iaftr)
      end if
C Check that the next row is adjacent to current row
C need to check additionally if next index is in the current row
      knext = kbeam(4, irnex)
      if ( kcurr + 1 .NE. knext .OR.
     ;     mod(kcurr, nrows) .EQ. 0) then
          ircur = irnex
          kcurr = kbeam(4, ircur)
          call rnext(kbeam, ircur, irnex, 0, nbeam)
          go to 801
      end if
C next row is adjecent, associate and move on
      call ascnx(kbeam, llasc, nbeam, 1, ircur, irnex, iaftr)

      ircur = irnex
      kcurr = kbeam(4, ircur)
      irnex = iaftr
      go to 801

C count the number of parents
  290 nclmp = 0
      do 2 i=1, nbeam
          lcurr = llasc(1, i)
          if ( lcurr .EQ. i ) then
              nclmp = nclmp + 1
              llasc(2, i) = nclmp
          end if
    2 continue
      end
