SUBROUTINE induced_vel_at_P

USE conpara
USE input
USE locat

IMPLICIT NONE

INTEGER :: i
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: rsarr,rarr, &
                                            dxarr,dyarr,dzarr
DOUBLE PRECISION:: r1s,r2s,rx1,rx2,dx1,dx2,dy1,dy2,dz1,dz2, &
                  rc,es,ls,rl,r1r2,r1xr2,r1yr2,r1zr2,rxrs, &
                  num,dem,dvix,dviy,dviz,gcm,vix,viy,viz
ALLOCATE(dxarr(nvvp),dyarr(nvvp),dzarr(nvvp), &
         rsarr(nvvp),rarr(nvvp))
!-------- Predetermine some parameters here ---------

rc = 0.001   !---Rc is read in as the source code shows
es = 0.000001

vix = 0.00   !
viy = 0.00   !---innitializes induced vel BEFORE loop...
viz = 0.00   !

!--------------- Determine R vectors ------------------------
DO i = 1,nvvp    
    dxarr(i) = xp - xm(i)
    dyarr(i) = yp - ym(i)
    dzarr(i) = zp - zm(i)

    r1s = dxarr(i)**2 + dyarr(i)**2 + dzarr(i)**2

    IF (r1s <= 0.0) THEN 
       rx1 = 1.0
       ELSE
       rx1 = SQRT(r1s)
    END IF 

    rsarr(i) = r1s 
    rarr(i) = rx1

END DO 

DO i = 1,nvvp-1

    dx1 = dxarr(i)
    dy1 = dyarr(i)
    dz1 = dzarr(i)
    r1s = rsarr(i)
    rx1 = rarr(i)

    dx2 = dxarr(i+1)
    dy2 = dyarr(i+1)
    dz2 = dzarr(i+1)
    r2s = rsarr(i+1)
    rx2 = rarr(i+1)

!!!!------ Get dot Product
    r1r2  = dx1*dx2 + dy1*dy2 + dz1*dz2
!!!!------ Get Cross Product Components
    r1xr2 = -1*(dy2*dz1 - dy1*dz2)
    r1yr2 = -1*(dx1*dz2 - dx2*dz1)
    r1zr2 = -1*(dx2*dy1 - dx1*dy2)

    ls   = r1s + r2s - 2.0*r1r2
    rl   = (rc**2)*(ls + es)
    rxrs = r1s*r2s - r1r2*r1r2
    
    num   =  (rx1 + rx2)*(rx1*rx2 - r1r2)
    dem   =  rx1*rx2*sqrt(rxrs**2 + rl**2)
    

!!!!--the 1 here (below) is =bld_inside() in the source
    gcm = gama*1*num/dem

    dvix = gcm*r1xr2
    dviy = gcm*r1yr2
    dviz = gcm*r1zr2

    vix  = vix + dvix
    viy  = viy + dviy
    viz  = viz + dviz
END DO 

WRITE(6,*) 'Biot-Savart induced velocity at point P:'
WRITE(6,*) '  Vix                       Viy                      Viz'
WRITE(6,*) vix,viy,viz

END SUBROUTINE induced_vel_at_P