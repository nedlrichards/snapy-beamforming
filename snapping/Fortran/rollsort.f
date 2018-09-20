      subroutine rollsort(ncurr, nrank, munrk, nchan, nsort)
      implicit none
C sort is efficent if ncurr is already almost sorted
C ncurr: nchan vector of current beamformer index values
C nrank: nchan vector of last ranking of beamformer index values
C nchan: number of channels in beamformer
C munrk: nsort sorted vector of channels which are no longer ranked
C nsort: number of unranked channels
      integer*4, intent(in)::nchan
      integer*4, intent(in)::ncurr(nchan), munrk(nsort)
      integer*4, intent(inout)::nsort, nrank(nchan)
      integer*4 :: i, j, jtemp
Cf2py intent(hide), depend(nchan) nrank
Cf2py intent(hide), depend(nchan) ncurr
Cf2py intent(hide), depend(nsort) nunrk

C roll unsorted channels to the end of the array

      do 1 i = 1, nsort
C if last channel in unsort, no need to roll
          if( munrk(nsort - i + 1).EQ.nchan) go to 1
C roll current unranked channel to the next avalible slot in back of rankings
          do 2 j = munrk(nsort - i + 1), (nchan - i)
             jtemp = nrank(j)
             nrank(j) = nrank(j + 1)
             nrank(j + 1) = jtemp
    2     continue
    1 continue

C work through unsorted channels and see if they need to be rolled back
C start from unsorted value closest to the front

C if array is totaly unsorted, there is no need to do anything to
C first index
      if( nsort.EQ.nchan ) nsort = nsort - 1
      do 3 i = nsort, 1, -1
C bring each value forward until they are sorted
          do 4 j = (nchan - i + 1), 2, -1
              if( ncurr(nrank(j - 1)).LE.ncurr(nrank(j)) ) then
                  go to 3
              else
                  jtemp = nrank(j - 1)
                  nrank(j - 1) = nrank(j)
                  nrank(j) = jtemp
              end if
    4     continue
    3 continue

      end
