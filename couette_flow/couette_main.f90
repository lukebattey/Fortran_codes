PROGRAM couette_main

USE variables_cf

    IMPLICIT NONE

OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) dt
READ(16,*) dy
READ(16,*) nu
READ(16,*) tau
READ(16,*) theta
READ(16,*) jmax
READ(16,*) nmax

ALLOCATE(Ynd(jmax),U(jmax),Unext(jmax),Uss(jmax))

!--------------- innitial values !!-----------------
U(1)    = 0.00
U(jmax) = 1.00 ! set as plate speed?
Ynd(1)    = 0.00
Ynd(jmax) = 1.00 ! ALWAYS 1...

Uss(1) = 0
Uss(jmax) = 1.00 ! Plate speed??

L = (jmax-1)*dy

DO j = 2,jmax-1
  Ynd(j)   = dy*(j-1)/L
  U(j) = Ynd(j) + SIN(pi*Ynd(j)) 
  Uss(j) = Uss(jmax)*(j-1)/(jmax-1)
END DO 

  tau = (L**2)/nu 
  dt0 = dt/tau
  dy0 = dy/L

  c = dt0/(dy0**2)  ! coef of "scheme"...

!---------------------- LOOPS! ------------------------------
IF (theta == 0) THEN 
!--------------------EXPLICIT SOLUTION!----------------------
!============================================================
  DO n = 1,nmax
    errs = 0

    DO j = 2,jmax-1
      Unext(j) = U(j) + c*(U(j+1)-2*U(j)+U(j-1))
      
      errs = (Unext(j) - Uss(j))**2
      serrs = serrs + errs
    END DO

    DO j = 2,jmax-1
      U(j) = Unext(j)
    END DO 

    RMSerr = SQRT(serrs/(jmax-2))
    
    IF (RMSerr <= 5e-7) THEN
      WRITE(6,*) n,'WOOOOOO'
      nmax = n
    END IF

  END DO
!----------------- COMBINED METHOD SOLUTION -------------------
ELSE
  WRITE(6,*) 'UNDER CONSTRUCTION. ABORTING....'
  EXIT
END IF
!--------------------WRITE RESULTS-----------------------------

DO j = 1,jmax
   WRITE(6,*) U(j),Ynd(j)
END DO



END PROGRAM couette_main