      subroutine asc2d(kbeam, llasc, nbeam, kclmp)
C Associate beams which are next to each other
      integer*4, intent(in) :: kbeam(4, nbeam)
      integer*4, intent(out) :: llasc(2, nbeam), kclmp


C preallocate association array
      do 1 i=1, nbeam
          llasc(1, i) = i
          llasc(2, i) = 0
    1 continue

C start loop at first index
      ircur = 1
      kcurr = kbeam(4, 1)

C find the first index of the next row
      call rnext(kbeam, ircur, irnex, 0, nbeam)

C 0 indicates we are at the end of the beams
  801 if ( irnex .LE. 0 ) go to 290

C Check that the next row is adjacent to current row
      knext = kbeam(4, irnex)
      if ( kcurr + 1 .NE. knext ) then
          ircur = irnex
          kcurr = kbeam(4, ircur)
          call rnext(kbeam, ircur, irnex, 0, nbeam)
          go to 801
      end if
C next row is adjecent, associate and move on
      call ascnx(kbeam, llasc, nbeam, 1, ircur, irnex, iafter)
      ircur = irnex
      kcurr = kbeam(4, ircur)
      irnex = iafter
      go to 801

C count the number of parents
  290 kclmp = 0
      do 2 i=1, nbeam
          lcurr = llasc(1, i)
          if ( lcurr .EQ. i ) then
              kclmp = kclmp + 1
              llasc(2, i) = kclmp
          end if
    2 continue
      end
