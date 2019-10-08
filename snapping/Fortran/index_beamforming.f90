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

subroutine cmerge(mshif, ichan, thresh, vbout, kbout, lcurr, lnmax, nchan, nnz)
  ! merge channels together and return values that are above a threshold
  ! this is a beamforming routine that takes advantage of the fact that
  ! the data is sorted channel by channel.
  !
  ! ::INPUTS::
  !
  ! mshif: nnz vector of value indices, already delay shifted
  ! ichan: nchan vector marking the start of each channels index and value vector
  ! nnz: total number of indecies
  ! nchan: number of data channels
  ! thresh: beamformer threshold. Only save values which sum GE this value
  ! lnmax: total allocated beamformer result positions

  integer, intent(in)  :: lnmax, nnz, nchan, thresh
  integer, intent(in)  :: mshif(nnz), ichan(nchan + 1)

  ! ::OUTPUTS::
  !
  ! vbout: vector of beamformer outputs GE thresh
  ! kbout: index of beamformer outputs

  integer, intent(out) :: kbout(lnmax), vbout(lnmax)

  ! ::COUNTERS::
  !

  integer :: lcurr

  ! Internal bookeeping variables to move through input index and value pairs
  ! lcurr: current position in beamformer result array
  ! mcurr: nchan vector of current index value on each channel
  ! mrank: nchan vector of ranking of the current index value on each channel
  ! mpos: nchan vector of position of current index relative to ichan

  integer :: mcurr(nchan), mrank(nchan), mpos(nchan)
  integer :: munrk(nchan), nsort

  ! setup for merge routine
  ! filling muncrk leads to a complete sort of mcurr
  ik = 1
  do i = 1, nchan
      mcurr(i) = mshif(ichan(i))
      mrank(i) = i
      munrk(i) = i
      mpos(i) = 0
  end do
  nsort = nchan

C sort the start of the buffer
      call rollsort(mcurr, mrank, munrk, nchan, nsort)
      jdb = 1

C merge step, combine vals at which indecies are equal
C initialize beamformer with value from first buffer
  801 nsort = 0
C beamformer accumulating variable
      tmp = vval(ichan(mrank(1)) + mpos(mrank(1)))
C delay of current beamformer result
      itmp = mcurr(mrank(1))
C step channel index
      mpre = mshif(ichan(mrank(1)) + mpos(mrank(1)))
      mpos(mrank(1)) = mpos(mrank(1)) + 1
C check against end condition
      if( ichan(mrank(1)) + mpos(mrank(1)) .GE.ichan(mrank(1) + 1) )
     ; then
          if( mpos(mrank(2)).GT.nnz ) then
              go to 299
C first channel is exhausted, but rest of channels are not yet
          else
              nsort = 1
              munrk(1) = 1
C nnz value of current index is used to mark an exhausted buffer
              mcurr(mrank(1)) = 2147483647
              mpos(mrank(1)) = nnz + 1
          end if
      else
C bookeep that a value was pulled
C check for a discontinuity in index values which requires a sort
        mcurr(mrank(1)) = mshif(ichan(mrank(1)) + mpos(mrank(1)))
        if( mpre - mcurr(mrank(1)) .LT. -1 ) then
            nsort = 1
            munrk(1) = 1
        end if

      end if

C add as many buffers as have equal index, then step each buffer
C this remains unchanged if there is no loop break in 3
      inum = nchan
      do 3 i = 2, nchan
C break loop if the current ranked channel is not the same as first channel
C nnz value of current index is used to mark an exhausted buffer
          if( mpos(mrank(i)).GT.nnz .OR. itmp.NE.mcurr(mrank(i)) ) then
              inum = i - 1
              goto 201
          end if
C accumulate value to beamformer
          tmp = tmp + vval(ichan(mrank(i)) + mpos(mrank(i)))
          mpre = mshif(ichan(mrank(i)) + mpos(mrank(i)))
          mpos(mrank(i)) = mpos(mrank(i)) + 1
C bookeep that a value was pulled
C check against end condition
          if( ichan(mrank(i)) + mpos(mrank(i)).LT.ichan(mrank(i) + 1) )
     ;    then
C bookeep that a value was pulled
              mcurr(mrank(i)) = mshif(ichan(mrank(i)) + mpos(mrank(i)))
C check for a discontinuity in index values which requires a sort
              if( mpre - mcurr(mrank(i)) .LT. -1 ) then
                  nsort = nsort + 1
                  munrk(nsort) = i
              end if
          else
C nnz value of current index is used to mark an exhausted buffer
              mcurr(mrank(i)) = 2147483647
              mpos(mrank(i)) = nnz + 1
              nsort = nsort + 1
              munrk(nsort) = i
          end if
    3 continue

C Check before adding value to beamformer result
  201 if( tmp.GE.thresh ) then
          if ( lcurr.EQ.lnmax ) goto 299
          lcurr = lcurr + 1
          vbout(lcurr) = tmp
          kbout(lcurr) = itmp
      end if

C Maintain the sorted order of the buffers
      if( nsort.GT.0 ) then
           if( jdb .LT. 4) then
C              print *, munrk
C              print *, nsort
           end if
           call rollsort(mcurr, mrank, munrk, nchan, nsort)
           if( jdb .LT. 4) then
C          print *, mcurr
C          print *, mrank
           end if
           jdb = jdb + 1
       end if

      go to 801

  299 end

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
