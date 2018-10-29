      subroutine ascnx(kbeam, llasc, nbeam, ioffs, ircur, irnex, iafter)
      implicit none
      integer*4, intent(in) :: nbeam, ioffs, ircur, irnex
      integer*4, intent(out) :: iafter
C Associate beams which are next to each other
      integer*4 :: kbeam(4, nbeam)
      integer*4 :: llasc(2, nbeam)

      integer*4 :: ltrgt, knext, icrtm, inxtm, ipar1, ipar2

C Establish channel index numbers
      knext = kbeam(4, irnex)
C Loop variables
      icrtm = ircur
      inxtm = irnex

C Loop untill current or next index changes
  801 if ( kbeam(4, icrtm) + ioffs .NE. kbeam(4, inxtm) ) go to 201
C Step the current row untill the last value is greater than next rows first
      if ( kbeam(2, icrtm) .LT. kbeam(1, inxtm) - 1 ) then
          icrtm = icrtm + 1
          go to 801
      end if
C Step the next row untill the last value is greater than next rows first
      if ( kbeam(2, inxtm) .LT. kbeam(1, icrtm) - 1 ) then
          inxtm = inxtm + 1
          go to 801
      end if
C Associate overlapping beams
      ipar1 = icrtm
      ltrgt = llasc(1, ipar1)
  816 if ( ltrgt .NE. ipar1 ) then
          ipar1 = ltrgt
          ltrgt = llasc(1, ipar1)
          go to 816
      end if

      ipar2 = inxtm
      ltrgt = llasc(1, ipar2)
  831 if ( ltrgt .NE. ipar2 ) then
          ipar2 = ltrgt
          ltrgt = llasc(1, ipar2)
          go to 831
      end if
C     llasc(1, ipar1) = llasc(1, ipar2)
      llasc(1, ipar2) = llasc(1, ipar1)

C Step channel with smallest max value
      if ( kbeam(2, icrtm) .LE. kbeam(2, inxtm) ) then
          icrtm = icrtm + 1
      else
          inxtm = inxtm + 1
      end if
      go to 801
C set iafter if it is not already

  201 if ( kbeam(4, icrtm) + ioffs - 1 .EQ. knext ) then
          call rnext(kbeam, inxtm, iafter, 0, nbeam)
      else
          iafter = inxtm
      end if

      end
