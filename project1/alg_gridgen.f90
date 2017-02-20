MODULE alg_gridgen
IMPLICIT NONE

    !------------------------------------------------------------
    !            DEFINE VARIABLES
    !------------------------------------------------------------

    INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(12)

    INTEGER,PARAMETER::imax=41,jmax=19,fi=11,di=31,ci=41
    INTEGER::i,j
    REAL(KIND=rDef),DIMENSION(1:imax,1:jmax)::x,y
    REAL(KIND=rDef),PARAMETER::fx=-.8, &
                    dx=1.00, &
                    cx=1.80, &
                    t=0.15, &
                    xint=1.008930411365

contains
    subroutine write_grid(a)

        INTEGER, INTENT(IN) :: a

        write(6,*) a

        !----------------- WRITE LOWER BOUNDARY (j=1) ---------------------

        DO i=1,fi-1
            x(i,1) = (-fx/(fi-1))*(i-1) + fx
            y(i,1) = 0
            write(6,*) i,x(i,1),y(i,1)
        END DO

        DO i=fi,di
            x(i,1) = (dx/(di-fi))*(i-fi)
            y(i,1) = 5*t*(0.2969*sqrt(xint*x(i,1)) &
                    -0.126*xint*x(i,1) &
                    -0.3516*((xint*x(i,1))**2) &
                    +0.2843*((xint*x(i,1))**3) &
                    -0.1015*((xint*x(i,1))**4)) 
            write(6,*) i,x(i,1),y(i,1)
        END DO

        DO i=di+1,ci
            x(i,1) = dx+((cx-dx)/(ci-di))*(i-di)
            y(i,1) = 0
            write(6,*) i,x(i,1),y(i,1)
        END DO

        !----------------- WRITE UPPER BOUNDARY (j=jmax) -------------------

        DO i=1,ci
            x(i,jmax) = fx+((cx-fx)/(imax-1))*(i-1)
            y(i,jmax) = 1.0
            WRITE(6,*) i,x(i,jmax),y(i,jmax)
        END DO

        print *, 3. * 6 / 11

    end subroutine write_grid
END MODULE alg_gridgen