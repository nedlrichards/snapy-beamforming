      subroutine smerge(mshif, vval, ichan, thresh, jchan, vbout,
     ;                  kbout, lcurr, lnmax, nchan, nnz)
      implicit none
C merge channels together and return values that are above a threshold
C this is a beamforming routine that takes advantage of the fact that
C the data is sorted channel by channel.
C
C ::INPUTS::
C
C mshif: nnz vector of value indices, already delay shifted
C vval: nnz vector of beamformer values
C ichan: nchan vector of start of each channels index and value vector
C thresh: beamformer threshold. Only save values which sum GE this value
C jchan: number of look vector, used to book keep across look directions
C nnz: total number of index and value pairs
C nchan: number of data channels
C lnmax: total allocated beamformer result positions
C
      integer*4, intent(in)  :: lnmax, nnz, nchan, jchan
      real*8, intent(in)     :: thresh
      integer*4, intent(in)  :: mshif(nnz), ichan(nchan + 1)
      real*8, intent(in)     :: vval(nnz)
C
C ::OUTPUTS::
C
C kbout: vector of beamformer indecies
C        (start index, end index, max index, channel number)
C vbout: vector of beamformer maxvalue
C
      integer*4, intent(out) :: kbout(4 * lnmax)
      real*8, intent(out)    :: vbout(lnmax)
C
C ::COUNTERS::
C
C lcurr: current position in beamformer result array
C
      integer*4, intent(inout) :: lcurr
C
C Internal bookeeping variables to move through input index and value pairs
C mcurr: nchan vector of current index value on each channel
C mrank: nchan vector of ranking of the current index value on each channel
C mpos: nchan vector of position of current index relative to ichan
C
      integer*4 :: mcurr(nchan), mrank(nchan), mpos(nchan)
      integer*4 :: munrk(nchan), nsort
      integer*4 :: mpre, itmp, inum, i
      integer*4 :: kend = 0, kfirst = 0, kmax = 0
      real*8    :: tmp
      real*8    :: vmax = 0
C
Cf2py intent(hide), depend(nnz) mshif
Cf2py intent(hide), depend(nnz) vval
Cf2py intent(hide), depend(nchan + 1) ichan
Cf2py intent(hide), depend(lnmax) kbout
Cf2py intent(hide), depend(lnmax) vbout
C
C
C setup for merge routine
C filling muncrk leads to a complete sort of mcurr
      do 1 i = 1, nchan
          mcurr(i) = mshif(ichan(i))
          mrank(i) = i
          munrk(i) = i
          mpos(i) = 0
    1 continue
      nsort = nchan

C sort the start of the buffer
      call rollsort(mcurr, mrank, munrk, nchan, nsort)

C kend of 0 means there is no beamformer output to save
      kend = 0
      kfirst = 0
      kmax = 0

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
              go to 298
C first channel is exhausted, but rest of channels are not yet
          else
              nsort = 1
              munrk(1) = 1
C overflow value of index is used to mark an exhausted buffer
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
          if ( itmp - kend .GT. 1 ) then
C Only save last beamformer result if the start index is not 0
              if ( kend.GT.0 ) then
              if (kend - kfirst .LT. 0) print *, kfirst
              if (kend - kfirst .LT. 0) print *, kend
                  kbout(lcurr * 4 - 3) = kfirst
                  kbout(lcurr * 4 - 2) = kend
                  kbout(lcurr * 4 - 1) = kmax
                  kbout(lcurr * 4) = jchan
                  vbout(lcurr) = vmax
              end if
C Increment beamformer result count, and initialize new cluster
              lcurr = lcurr + 1
              kfirst = itmp
              kend = itmp
              kmax = itmp
              vmax = tmp
          else
              kend = itmp
              if ( tmp.GT.vmax ) then
                  vmax = tmp
                  kmax = itmp
              end if
          end if
      end if

C Maintain the sorted order of the buffers
      if( nsort.GT.0 ) then
          call rollsort(mcurr, mrank, munrk, nchan, nsort)
      end if

      go to 801

C check if there is one more buffer that needs to be saved
  298 if ( kend.GT.0 ) then
          kbout(lcurr * 4 - 3) = kfirst
          kbout(lcurr * 4 - 2) = kend
          kbout(lcurr * 4 - 1) = kmax
          kbout(lcurr * 4) = jchan
          vbout(lcurr) = vmax
      end if

  299 end
