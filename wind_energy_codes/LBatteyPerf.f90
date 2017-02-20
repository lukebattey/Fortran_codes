PROGRAM LBatteyPerf
    IMPLICIT NONE

INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(10)
REAL(KIND=rDef),PARAMETER :: pi = 4.0*ATAN(1.0)
INTEGER :: vmax,vmin,vstep,V,B,i,af,imax,numaf,nAl,n,nmax
CHARACTER(len=90),ALLOCATABLE,DIMENSION(:) :: airfile 
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:) :: rR,cR,tw,Cnor,Ctan,P
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:) :: alpha,cl,cd
INTEGER,ALLOCATABLE,DIMENSION(:) :: ai,nAlM 
REAL(KIND=rDef):: R,rho,rpm,pitch,numjunk,omega,wSec,rtcut,lam,sig, &
                  F,A1,B1,phir,Lp,Dp,Psec,Vn,Vt,phi,beta,iAl,cli,cdi, &
                  a,aLast,ap,apLast

OPEN(16,FILE = 'lbatteyPerf.inp', FORM = 'FORMATTED')
READ(16,*) R
READ(16,*) rtcut
READ(16,*) B
READ(16,*) rho
READ(16,*) pitch
READ(16,*) rpm
READ(16,*) vmin,vmax,vstep
ALLOCATE(P(vmax-vmin+2))
READ(16,*) numaf
ALLOCATE(airfile(numaf),nAlM(numaf))

DO af=1,numaf
    READ(16,*) airfile(af)
END DO
READ(16,*) imax 
ALLOCATE(rR(imax),cR(imax),tw(imax),ai(imax),Cnor(imax),Ctan(imax))

DO i=1,imax
    READ(16,*) rR(i),cR(i),tw(i),ai(i)
END DO

DO af=1,numaf
    OPEN(af*10,FILE = airfile(af), FORM = 'FORMATTED')
    READ(af*10,*) nAlM(af)
    ALLOCATE(alpha(nAlM(af),numaf),cl(nAlM(af),numaf),cd(nAlM(af),numaf))
    DO nAl = 1,nAlM(af)
        READ(af*10,*) alpha(nAl,af),cl(nAl,af),cd(nAl,af) 
    END DO
END DO
omega = rpm*2.0*pi/60.0
wSec = (R-rtcut*R)/imax
nmax = 10      
!-------------------------------------------------------------------------
!===================== MAIN LOOP STARTS HERE =============================
!-------------------------------------------------------------------------
DO V = vmin,vmax,vstep
    DO i = 1,imax
        !--------------------- Initialize a and ap -----------------------
        lam = omega*rR(i)*R/V
        sig = (B*cR(i)*R)/(2*pi*rR(i)*R)
        beta = tw(i) + pitch
        aLast = 0.25*(2.0+(pi*lam*sig)-SQRT(4-(4*pi*lam*sig) &
               +(pi*(lam**2)*sig)*((8*pi*beta/180.0)+(pi*sig))))
        apLast = 0
        !-------------------- iteratively solve for a and ap -------------
        DO n = 1,nmax
            Vn = V*(1-aLast)
            Vt = omega*rR(i)*R*(1+apLast)
            phi = 180*(ATAN2(Vn,Vt))/pi
            beta = tw(i) + pitch
            iAl = phi - beta
            !--------- interpolate airfoil file for cl and cd ------------
            af = ai(i)
            DO nAl = 1,nAlM(af)-1
                IF (alpha(nAl+1,af) > iAl .AND. alpha(nAl,af) < iAl) THEN

                  cli = cl(nAl,af)+((cl(nAl+1,af)-cl(nAl,af))/ &
                    (alpha(nAl+1,af)-alpha(nAl,af)))*(iAl-alpha(nAl,af))

                  cdi = cd(nAl,af)+((cd(nAl+1,af)-cd(nAl,af))/ &
                    (alpha(nAl+1,af)-alpha(nAl,af)))*(iAl-alpha(nAl,af))
                  EXIT
                END IF
            END DO
            !---------- Redetermine a and ap (with cl and cd) ------------
            phir = pi*phi/180 
            F = 2.0*ACOS(EXP(-(B*R*(1-rR(i)))/(2*R*rR(i)*SIN(phir))))/pi
            A1 = 4.0*F*((SIN(phir))**2)
            B1 = sig*(cli*COS(phir)+cdi*SIN(phir))

            a = (1.0+A1/b1)**(-1.0)
            ap = (-1.0+(4.0*F *SIN(phir)*COS(phir))/ &
            (sig*(cli*SIN(phir)-cdi*COS(phir))))**(-1.0)

            IF( (ABS(a-aLast)+ABS(ap-apLast)) <= 1.0e-9) THEN 
                EXIT
            END IF
            aLast = a
            apLast = ap
        END DO
        !------------- Find Coefs of lift and Total Power ----------------
        Lp = 0.5*rho*(Vn**2 + Vt**2)*cR(i)*R*cli
        Dp = 0.5*rho*(Vn**2 + Vt**2)*cR(i)*R*cdi
        
        Psec = (Lp*SIN(Phir)-Dp*COS(Phir))*rR(i)*R*omega*wSec
        P(V) = P(V) + B*Psec

        IF (V == 7) THEN
            Ctan(i)=cli*SIN(pi*iAl/180)-cdi*COS(pi*iAl/180)
            Cnor(i)=cli*COS(pi*iAl/180)+cdi*SIN(pi*iAl/180)
        END IF
    END DO
END DO  
!=========================================================================
!------------------------- END OF MAIN LOOP ------------------------------
!=========================================================================

WRITE(6,*) 'Total Power Output at Wind Speeds:'
WRITE(6,*) 'V (m/s)     P (kW)'
DO v = vmin,vmax,vstep
    WRITE(6,'(I5,F15.5)') V,P(V)/1000.0
END DO

WRITE(6,*) 'Sectonal Force Coefs at r/R Values for V = 7 m/s'
WRITE(6,*) 'r/R     C_normal     C_tangential'
DO i = 1,imax
    WRITE(6,'(F5.3,F12.5,F13.5)') rR(i),Cnor(i),Ctan(i)
END DO
END PROGRAM LBatteyPerf