subroutine adddetection(nchan, detection, norm, alldetections, alli, allnorms, numdetect)
! add detection to detection list
    implicit none
    integer :: nchan, i, j, norm, numdetect, normcomp, detecti
    integer :: alldetections(nchan, numdetect), alli(numdetect), allnorms(numdetect)

    do i = 1, numdetect
      normcomp = allnorms(i)
      detecti = alli(1)
      if (detecti == 0) exit




end subroutine

subroutine numsame(nchan, detection1, detection2, numsame)
! Compute number of identical channels between 2 detections
    implicit none
    integer :: nchan, i, numsame
    integer :: detection1(nchan), detection2(nchan)

    numsame = 0
    do i = 1, nchan
      if (detection1(i) == detection2(i)) numsame = numsame + 1
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
