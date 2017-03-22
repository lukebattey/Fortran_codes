PROGRAM s_sonic_main

USE variables_ss

    IMPLICIT NONE

!====== READ IN GRID AND GET METRICS / JACOBIANS ==========
CALL grid_read_metrics

!====== SET INITIAL CONDITIONS AND BCs ====================
CALL init_L_BC   

!====== 




END PROGRAM s_sonic_main