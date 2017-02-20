PROGRAM main

USE input
USE locat
USE conpara

IMPLICIT NONE

CALL readin
CALL induced_vel_at_P
WRITE(*,*)'done!'

END PROGRAM main