      subroutine sbeam(ival, vval, ichan, thresh, btaus,
     ;                 kbout, vbout, nbeam,
     ;                 nbmax, nval, nchan, nlook)
      implicit none
C
C Beamformer for sparse input data
C
C Loops through a delay hypothesis matrix and beamform for each vector
C The result of the beamformer is expected to be sparse, so some
C generous guess at the return size is made.
C
C ::INPUTS::
C
C thresh: beamformer value required before saving result
C nbmax: a priori upper bound of number of beamformer excedences
C nval: number of data points
C nchan: number of channels of data
C nlook: number of beamformer look delay vectors
C vval: flat vector of data values
C ival: flat vector of time index of data
C ichan: position of channel transistions in ival and vvavl
C btaus: nchan x nlook matrix of beamformer delays rounded to integrer value
C
      real*8, intent(in)      :: thresh
      integer*4, intent(in)   :: nbmax, nval, nchan, nlook
      real*8, intent(in)      :: vval(nval)
      integer*4, intent(in)   :: ival(nval), ichan(nchan)
      integer*4, intent(in)   :: btaus(nchan, nlook)
C
C ::OUTPUTS::
C
C nbeam: a posteriori number of beamformer outputs
C vbout: value of each beamformer output
C kbout: 4 x nbmax matrix of beamformer indicies
C        (start index, end index, max index, channel number)
C
      integer*4, intent(out)  :: kbout(4, nbmax), nbeam
      real*8, intent(out)     :: vbout(nbmax)
C
C ::INTERAL::
C
      integer*4               :: mshif(nval), llast, jchan(nchan + 1)
      integer*4               :: i
C
Cf2py intent(hide), depend(nval) vval
Cf2py intent(hide), depend(nval) ival
Cf2py intent(hide), depend(nchan) ichan
C

C Loop over outside index of beamformer matrix, accumulate results in
C vbout and kbout vectors
      nbeam = 0
C jchan is i chan, but with a last index of nval + 1
      do 9 i = 1, nchan
          jchan(i) = ichan(i)
    9 continue
      jchan(nchan + 1) = nval + 1

      do 1 i = 1, nlook
          llast = nbeam
          call addtau(jchan, btaus(1, i), ival, mshif, nchan, nval)
          call smerge(mshif, vval, jchan, thresh, i,
     ;                vbout, kbout, nbeam, nbmax, nchan, nval)
    1 continue

      end
