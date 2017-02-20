SUBROUTINE readin
! This subroutine does the following
! a) Read Input File

USE input
USE locat

IMPLICIT NONE

INTEGER :: i
NAMELIST /rotors/ nbld,nseg,segprev,tolnseg
NAMELIST /conditions/ adratio,ct,gama
NAMELIST /location/ xp,yp,zp,nvvp

OPEN(UNIT = 5, FILE = 'input.inp')
READ(5,rotors)        
READ(5,conditions)    
READ(5,location)
CLOSE(5)

ALLOCATE(xm(nvvp),ym(nvvp),zm(nvvp))
OPEN(UNIT = 10, FILE = 'locations.inp')

READ(10,*)  ! Read the first line of locations.inp
DO i = 1,nvvp     
  READ(10,*) xm(i),ym(i),zm(i)
END DO 
CLOSE(10)

END SUBROUTINE readin