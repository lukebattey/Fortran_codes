MODULE conpara        ! Constant paparmeter
REAL,PARAMETER :: pi = 4.0*ATAN(1.0)
END MODULE conpara

MODULE input          ! Input file paparmeter
INTEGER nbld, nseg, segprev, tolnseg, nvvp
DOUBLE PRECISION adratio, ct, gama, xp, yp, zp
END MODULE input

MODULE locat          ! Input location 
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::xm,ym,zm
END MODULE locat