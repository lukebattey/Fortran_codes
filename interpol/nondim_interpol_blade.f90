PROGRAM nondim_interpol_blade
	IMPLICIT NONE

!------------------------------------------------------------
!			DEFINE VARIABLES
!------------------------------------------------------------

INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(12)
INTEGER,PARAMETER:: imax=16
INTEGER::i,numrR,il
REAL(KIND=rDef),DIMENSION(1:imax) :: r, rR, rRi, c, cR, cRi, tw, twi
REAL(KIND=rDef),PARAMETER:: Rmax = 5.03, &
							rRmin = 0.275, &
							rRmax = 0.975
REAL(KIND=rDef):: rRinc
numrR = 15

rRinc = (rRmax - rRmin)/(numrR - 1)


OPEN(16,FILE = 'Blade_geo_tab.txt', FORM = 'FORMATTED')
OPEN(26,FILE = 'Nondim_blade.txt', FORM = 'FORMATTED')
OPEN(36,FILE = '05_step_blade_geo.txt', FORM = 'FORMATTED')

DO i=1,imax
	READ(16,*) r(i),c(i),tw(i)
	rR(i) = r(i)/Rmax
	cR(i) = c(i)/Rmax
	WRITE(26,*) rR(i),cR(i),tw(i)
END DO

!--- CHECK IF INTERPOLATION IS POSSIBLE --------

IF (rR(1)<rRmin .AND. rR(imax)>rRmax) THEN
ELSE
WRITE(6,*) 'ERROR! WEE-WOO WEE-WOO (siren sound) PROGRAM TERMINATED!'
STOP
END IF

!--- INTERPOLATION! with nested loop :( ---------

DO i=1,numrR
    rRi(i) = rRmin + rRinc*(i*1.0-1.0)
    DO il=1,imax-1
        IF (rR(il) < rRi(i) .AND. rR(il+1) > rRi(i)) THEN
            cRi(i) = cR(il) + ((cR(il+1)-cR(il))/(rR(il+1)-rR(il)))*(rRi(i)-rR(il))
            twi(i) = tw(il) + ((tw(il+1)-tw(il))/(rR(il+1)-rR(il)))*(rRi(i)-rR(il))
        END IF
    END DO
    WRITE(36,*) rRi(i),cRi(i),twi(i)
END DO

END PROGRAM nondim_interpol_blade