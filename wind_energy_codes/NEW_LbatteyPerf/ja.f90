program ja
    implicit none

    INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(10)
    type realArr
        real(kind=rDef),allocatable :: el(:)
    endtype realArr
    type(realArr),allocatable,dimension(:) :: alpha, cl, cd
    integer::numaf, numang, af, nAl

    ! read # airfoils from file e.g. 2
    numaf = 2

    allocate(alpha(numaf), cl(numaf), cd(numaf))

    do af = 1, numaf
        ! read num angles from file
        numang = 57

        allocate(alpha(af)%el(numang), cl(af)%el(numang), cd(af)%el(numang))
        do nAl = 1, numang
            ! read from file...
            alpha(af)%el(nAl) = 7.7777
            cl(af)%el(nAl) = 5.5555
        end do

    end do

end program ja