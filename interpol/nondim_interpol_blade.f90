PROGRAM nondim_interpol_blade
	IMPLICIT NONE

!------------------------------------------------------------
!			DEFINE VARIABLES
!------------------------------------------------------------

INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(12)

INTEGER,PARAMETER:: imax=16
INTEGER::i
REAL(KIND=rDef),DIMENSION(1:imax):: r, rR, rRi, c, cR, cRi, tw, twi
REAL(KIND=rDef),PARAMETER:: Rmax = 5.03

OPEN(16,FILE = 'Blade_geo_tab.txt', FORM = 'FORMATTED')
OPEN(26,FILE = 'Nondim_blade.txt', FORM = 'FORMATTED')

DO i=1,imax
	READ(16,*) r(i),c(i),tw(i)
	rR(i) = r(i)/Rmax
	cR(i) = c(i)/Rmax
	WRITE(26,*) rR(i),cR(i)
END DO

iLower = 1

! make sure rR(1) is less than your first radial location (0.275)
! make sure rR(iMax) is greater than your last radial location (1)
! if not, print error and quit

DO i=1,15
	rRi(i) = 0.275+0.05*(i-1)
	! Increment iLower until rR(iLower)
	DO WHILE (rRi(i) < rR(iLower))
		iLower = iLower + 1
	END DO
	! Interpolate between rR(iLower) and rR(iLower + 1)
!	WRITE(6,*) rR(i),cR(i)
	WRITE(6,*) rRi(i)
	IF (rR(i) <= rRi(i) .AND. rR(i+1) >= rRi(i)) THEN
		cRi(i) = cR(i) + ((cR(i+1)-cR(i))/(rR(i+1)-rR(i)))*(rRi(i)-rR(i))
	
	!ELSE IF ()

	ELSE 
!	WRITE(6,*) "FUUUUUUUCK"

	END IF
	!WRITE(6,*) rRi(i)!,cRi(i)
END DO







END PROGRAM nondim_interpol_blade