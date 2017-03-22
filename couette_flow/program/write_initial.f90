PROGRAM write_init_couette

USE variables_cf

    IMPLICIT NONE

OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
READ(16,*) dt0
READ(16,*) jmax
READ(16,*) theta

ALLOCATE(Ynd(jmax),U(jmax))

OPEN(26,FILE = infile, FORM = 'FORMATTED')
WRITE(26,*) jmax

DO j = 1,jmax
    Ynd(j)   = 1.0*(j-1)/(jmax-1)
    U(j) = Ynd(j) + SIN(pi*Ynd(j))  !<-- You could make this
                                    !    any function! (in reality)
    WRITE(26,*) U(j),Ynd(j)                
END DO

END PROGRAM write_init_couette
