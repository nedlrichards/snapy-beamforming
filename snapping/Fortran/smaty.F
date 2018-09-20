#include "fintrf.h"
C
C This is an example of the FORTRAN code required for interfacing
C a .MEX file to MATLAB.
C
C This subroutine is the main gateway to MATLAB.  When a MEX function
C  is executed MATLAB calls the MEXFUNCTION subroutine in the corresponding
C  MEX file.  
C
C Copyright 1984-2006 The MathWorks, Inc.
C 
C

      SUBROUTINE MEXFUNCTION(NLHS, PLHS, NRHS, PRHS)
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER PLHS(*), PRHS(*)

C-----------------------------------------------------------------------
C

      INTEGER NLHS, NRHS
C
C-----------------------------------------------------------------------
C     (pointer) Replace integer by integer*8 on 64-bit platforms
C
      MWPOINTER MXCREATEDOUBLEMATRIX, MXGETPR
      MWPOINTER mexprintf, mxCreateDoubleScalar
      mwPointer mxCreateNumericArray
      integer*4 mxClassIDFromClassName
C
C-----------------------------------------------------------------------
C


      MWSIZE MXGETM, MXGETN
      mwSize dims(2),mxrow, mxcol, mxval, mxcha, mxlab, mx1, mx2, mx3
      MWPOINTER POUT1, POUT2, POUT3, N6
C
C KEEP THE ABOVE SUBROUTINE, ARGUMENT, AND FUNCTION DECLARATIONS FOR USE
C IN ALL YOUR FORTRAN MEX FILES.

C     Arguments for computational routine:
      real*8, allocatable :: vvals(:, :), vclmp(:, :)
      real*8 thresh
      integer*4, allocatable :: ivals(:, :), ichan(:, :)
      integer*4, allocatable :: btaus(:, :), lbout(:, :)
      integer*4 nbmax, nlbin, c0
      integer*4 nrows, ncolm, nchan, nvals, nlabl

C
C CHECK FOR PROPER NUMBER OF ARGUMENTS
C
      IF (NRHS .NE. 9) THEN
        CALL MEXERRMSGTXT('SMATY requires nine input arguments')
      ELSEIF (NLHS .NE. 3) THEN
        CALL MEXERRMSGTXT('SMATY requires three output arguments')
      ENDIF

C     Get the size of the input array. create both mx and int4 variables
      mx1 = 1
      mx2 = 2
      mx3 = 3
      c0 = 0
      mxval = mxGetM(prhs(1))
      nvals = mxval
      mxcha = mxGetM(prhs(3))
      nchan = mxcha

      allocate(ivals(mxval, mx1))
      allocate(vvals(mxval, mx1))
      allocate(ichan(mxcha, mx1))

      call mxCopyPtrToInteger4(mxGetPr(prhs(8)), nrows, mx1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(9)), ncolm, mx1)
      mxrow = nrows
      mxcol = ncolm
      allocate(btaus(mxcha, mxrow * mxcol))

      call mxCopyPtrToInteger4(mxGetPr(prhs(7)), nlbin, mx1)
      mxlab = nlbin

C     Create matrix for the return argument.
      dims(1) = mx3
      dims(2) = mxlab
      plhs(1) = mxCreateNumericArray(mx2, dims, 
     +          mxClassIDFromClassName('int32'), 
     +          c0)
      plhs(2) = mxCreateDoubleMatrix(mxlab, mx1, c0)
      dims(1) = mx1
      dims(2) = mx1
      plhs(3) = mxCreateNumericArray(mx2, dims, 
     +          mxClassIDFromClassName('int32'), 
     +          c0)
C     Create Fortran array from the input argument.


      call mxCopyPtrToInteger4(mxGetPr(prhs(1)), ivals, mxval)
      call mxCopyPtrToReal8(mxGetPr(prhs(2)), vvals, mxval)
      call mxCopyPtrToInteger4(mxGetPr(prhs(3)), ichan, mxcha)
      call mxCopyPtrToReal8(mxGetPr(prhs(4)), thresh, mx1)
      call mxCopyPtrToInteger4(mxGetPr(prhs(5)), btaus, 
     ;                         mxrow * mxcol * mxcha)
      call mxCopyPtrToInteger4(mxGetPr(prhs(6)), nbmax, mx1)

      allocate(vclmp(mx1, mxlab))
      allocate(lbout(mx3, mxlab))

C     Create pointers
      POUT1 = mxGetPr(plhs(1))
      POUT2 = mxGetPr(plhs(2))
      POUT3 = mxGetPr(plhs(3))

      call mxCopyPtrToInteger4(POUT1, lbout, mx3 * mxlab)
      call mxCopyPtrToReal8(POUT2, vclmp, mxlab)
      call mxCopyPtrToInteger4(POUT3, nlabl, mx1)

      call snapy(ivals, vvals, ichan, thresh, btaus, nbmax, lbout,
     ;           vclmp, nlabl, nlbin, nvals, nchan, ncolm, nrows)
      if (nlabl .EQ. 0) then
      call mexErrMsgIdAndTxt('MATLAB:revord:nInput',
     ;                        'Allocate more beams')
      return
      else if (nlabl .GE. nlbin) then
      call mexErrMsgIdAndTxt('MATLAB:revord:nInput',
     ;                       'Allocate more labels')
      return
      end if

C     Load the data into y_ptr, which is the output to MATLAB.
      call mxCopyInteger4ToPtr(lbout, POUT1, mx3 * mxlab)
      call mxCopyReal8ToPtr(vclmp, POUT2, mxlab)
      call mxCopyInteger4ToPtr(nlabl, POUT3, mx1)

      return
      end