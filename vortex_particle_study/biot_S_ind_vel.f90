PROGRAM induced_vel_at_P
    IMPLICIT NONE

INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(10)
REAL(KIND=rDef),PARAMETER :: pi=4.0*ATAN(1.0),root2=2**0.5
INTEGER :: nfil,j,NumPpv
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:)::xm,ym,zm, &
                                          dxarr,dyarr,dzarr, &
                                          rsarr,rarr
REAL(KIND=rDef):: gama,xp,yp,zp,r1s,r2s,rx1,rx2, &
                  dx1,dx2,dy1,dy2,dz1,dz2, &
                  rc,es,ls,rl,r1r2,r1xr2,r1yr2,r1zr2,rxrs, &
                  num,dem,dvix,dviy,dviz,gcm,vix,viy,viz,sigj, & 
                  Csig,hRes

OPEN(16,FILE = 'vector_and_point.inp', FORM = 'FORMATTED')

READ(16,*) gama
READ(16,*) Csig
READ(16,*) hRes
READ(16,*) xp,yp,zp
READ(16,*) NumPpv
READ(16,*) nfil
ALLOCATE(xm(nfil),ym(nfil),zm(nfil), &
         dxarr(nfil),dyarr(nfil),dzarr(nfil), &
         rsarr(nfil),rarr(nfil))

DO j = 1,nfil
     READ(16,*) xm(j),ym(j),zm(j)
END DO 

!-------- Predetermine some parameters here ---------

rc = 0.001   !--- Rc is read in as the source code shows
es = 0.000001

vix = 0.00   !
viy = 0.00   !--- innitializes induced vel BEFORE loop...
viz = 0.00   !

!--------------- Determine R vectors ------------------------
DO j = 1,nfil
    dxarr(j) = xp - xm(j)
    dyarr(j) = yp - ym(j)
    dzarr(j) = zp - zm(j)

    r1s = dxarr(j)**2 + dyarr(j)**2 + dzarr(j)**2

    IF (r1s <= 0.0) THEN 
       rx1 = 1.0
       ELSE
       rx1 = SQRT(r1s)
    END IF 

    rsarr(j) = r1s 
    rarr(j) = rx1

END DO 

DO j = 1,nfil-1

    dx1 = dxarr(j)
    dy1 = dyarr(j)
    dz1 = dzarr(j)
    r1s = rsarr(j)
    rx1 = rarr(j)

    dx2 = dxarr(j+1)
    dy2 = dyarr(j+1)
    dz2 = dzarr(j+1)
    r2s = rsarr(j+1)
    rx2 = rarr(j+1)

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
!!!!------ Also (gama here is "gamu" in sounce, which comes pre-divided by 4*pi)
    gcm = gama*1*num/(dem*4*pi)

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


END PROGRAM induced_vel_at_P