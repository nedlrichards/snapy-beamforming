module index_beamforming

use types, only: dp
use utils, only: stop_error

implicit none

subroutine adddetection(nchan, maxdetect, countweight, numdetect, detection, &
                        norm, chancount, detectL, idL, normL, countL, exception)
  ! add detection to detection list
  ! arguments
  implicit none
  integer :: nchan, maxdetect, numdetect
  real(dp) :: countweight
  integer :: detection(nchan), norm, chancount
  integer :: detectL(nchan, maxdetect), idL(maxdetect)
  integer :: normL(maxdetect), countL(maxdetect)
  logical :: exception
  ! local
  real(dp) :: score_old, score_new
  integer :: normcomp, countcomp
  logical :: ismatch = .false.
  logical :: isnew = .false.

  ! look for a matching detection in list, if not determine where it goes
  finder: do i = 1, numdetect
    detecti = idL(i)
  end do

  ! Decide if it is necassary to add detection to list
  decider: if ismatch then
    score_new = norm + chancount * countcomp
    normcomp = normL(i)
    countcomp = countL(i)
    score_old = normcomp + countweight * countcomp
    if (score_new .GT. score_old) isnew = .true.
  end if


  ! Add detection to list if needed and if there is room
  adder: if isnew then
    ! make sure there is room for a new detection
    if (numdetect == maxdetect) then
      exception = .true.
      exit outer
    else
      numdetect = numdetect + 1
    end if
    ! base case, append detection to end of list
    if (i == numdetect) then
      detectL(numdetect) = detection
      idL(numdetect) = detection(1)
      normL(numdetect) = norm
      countL(numdetect) = chancount
    ! insert detection into list
    else
      ! move all list elements to create room for detection
      do j = detecti + 1, numdetect
        detectL(j) = detectL(j - 1)
        idL(j) = idL(1, j - 1)
        normL(j) = normL(j - 1)
        chancount(j) = chancount(j - 1)
      end do
      ! insert new detection
      detectL(detecti) = detection
      idL(detecti) = detection(1)
      normL(detecti) = norm
      countL(detecti) = chancount
    end if
  end if
end subroutine

subroutine numsame(nchan, detection1, detection2, samecount)
  ! Compute number of identical channels between 2 detections
  implicit none
  integer :: nchan, i, samec
  integer :: detection1(nchan), detection2(nchan)

  samecount = 0
  do i = 1, nchan
    if (detection1(i) == detection2(i)) samecount = samecount + 1
  end do
end subroutine

subroutine l1norm(nchan, arrs, chani, l1)
  ! compute l1 norm of arrival time difference from median time
  implicit none
  integer :: nchan, i, medval, l1
  integer :: arrs(nchan), chani(nchan)

  medval = arrs(chani(nchan / 2))
  l1 = 0

  do i = 1, nchan
    l1 = l1 + abs(arrs(i) - medval)
  end do
end subroutine

subroutine isort(nchan, arrs, chani)
  !   insertion sort used to determine next row to pop detection from
  implicit none
  integer :: nchan, i, j, x, cc
  integer :: arrs(nchan), chani(nchan)

  do i = 2, nchan
    x = arrs(i)
    cc = chani(i)
    j = i - 1
    do while (j >= 1)
      if (arrs(j) <= x) exit
      arrs(j + 1) = arrs(j)
      chani(j + 1) = chani(j)
      j = j - 1
    end do
    arrs(j + 1) = x
    chani(j + 1) = cc
  end do
end subroutine

end module
