program type_example
    implicit none

    type point
        real x
        real y
    end type

    type(point) some_point
    type(point) another_point

    some_point%x = 5.0
    some_point%y = 10.0

    write(6, *) some_point%x, some_point%y

    another_point = some_point

    write(6, *) another_point%x, another_point%y

end program type_example