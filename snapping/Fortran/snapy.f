      subroutine snapy(ival, vval, ichan, thresh, btaus, nbmax, labls,
     ;                 vdetc, nlabl, nlbin, nval, nchan, nrows, ncolm)
      implicit none
C
C Detection and localization for impulsive sources
C
C Multi-channel data is assumed to be sparse, and efficently described as
C (time index, data value) pairs. The data is assumed sorted by time index for
C each channel. Possible source locations are searched using time-domain
C beamforming, creating a matrix with dimension of time and location.
C
C This two-dimensional matrix is then clustered by separation in time and
C location, and the location clustering can be one or two dimensional. One-
C dimensional clustering is used if the beamforming look vectors represent a
C one-dimensional search, such as source angle. Two-dimensional clustering is
C used if the beamforming look vectors represent a two-dimensional search, such
C as (source angle, source range). The clustering result is also sparse,
C recording source location, value and number of detections in cluster.
C
C ::INPUTS::
C
C thresh: beamformer value required before saving result
C nbmax: a priori upper bound of number of beamformer excedences
C nlbin: a priori upper bound of number of dections (labels)
C nval: number of data points
C nchan: number of channels of data
C nrows: number of beamformer looks in dimension 1
C ncols: number of beamformer looks in dimension 2
C vval: flat vector of data values
C ival: flat vector of time index of data
C ichan: position of channel transistions in ival and vvavl
C btaus: nchan x nlook matrix of beamformer delays rounded to integrer value
C
      real*8, intent(in)      :: thresh
      integer*4, intent(in)   :: nbmax, nlbin, nval, nchan, nrows, ncolm
      real*8, intent(in)      :: vval(nval)
      integer*4, intent(in)   :: ival(nval), ichan(nchan)
      integer*4, intent(in)   :: btaus(nchan, nrows * ncolm)
C
C ::OUTPUTS::
C
C nlabl: number of dections (labels)
C labls: index location of dections (labels)
C        (index dimension 1, index dimension 2, number of beamformer dections)
C vdetc: detection value, equal to max beamformer output in dection cluster
C
      integer*4, intent(out)  :: nlabl, labls(3, nlbin)
      real*8, intent(out)     :: vdetc(nlbin)
C
C ::INTERNAL
C
      real*8                  :: vbeam(nbmax)
      integer*4               :: kbeam(4, nbmax), nbeam, nlook,
     ;                           llasc(3, nbmax)

      nlabl=0
      nlook = nrows * ncolm
C Sparse data beamforming
      call sbeam(ival, vval, ichan, thresh, btaus, kbeam, vbeam, nbeam,
     ;           nbmax, nval, nchan, nlook)
      if ( nbeam.EQ.nbmax ) then
          nlabl = -1
          return
      end if
C three-dimsional beamformer output labeling
      call asc3d(kbeam, llasc, nrows, nbeam, nlabl)
C Accumulate labeled clusters into sparse detection outputs
      if ( nlabl .GE. nlbin) return
      call clump(kbeam, vbeam, llasc, labls, vdetc, nbeam, nlabl)

      end
