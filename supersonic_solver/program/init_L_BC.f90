SUBROUTINE init_and_LR_BCs

USE variables_ss

IMPLICIT NONE

OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
READ(16,*) ui
READ(16,*) rhoi
READ(16,*) gama
CLOSE(16) 

WRITE(6,*) ui,rhoi,gama

ALLOCATE(u(imax,jmax),v(imax,jmax),p(imax,jmax),rho(imax,jmax))

ALLOCATE(Ust(imax,jmax,4)) ! 4 to account for each component
                           ! of the state-vector (Ust)

END SUBROUTINE init_and_LR_BCs