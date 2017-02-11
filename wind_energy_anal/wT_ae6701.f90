PROGRAM wT_ae6701
    IMPLICIT NONE

INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(12)
REAL(KIND=rDef),PARAMETER :: pi = 4.0*ATAN(1.0)
REAL(KIND=rDef):: R,rho,rpm,pitch
INTEGER :: vmax,vmin,vstep,V,B
INTEGER :: i,af,imax,numaf,nAlpha,n,nmax
CHARACTER(len=90),ALLOCATABLE,DIMENSION(:) :: airfile 
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:) :: rR,cR,tw,Vn,Vt,phi,beta,iAlpha,cli,cdi, &
                                            Lp,Dp,Psec,a,ap,phir,Cnor,Ctan,P
REAL,ALLOCATABLE,DIMENSION(:,:) :: alpha,cl,cd
INTEGER,ALLOCATABLE,DIMENSION(:) :: ai,nAlphaMax 
CHARACTER(len=1) :: junk
REAL :: numjunk,omega,wSec,rootcut,lam,sig,F,A1,B1

OPEN(16,FILE = 'ae6701_WT.inp', FORM = 'FORMATTED')

READ(16,*) R
READ(16,*) rootcut
READ(16,*) B
READ(16,*) rho
READ(16,*) pitch
READ(16,*) rpm
READ(16,*) vmin,vmax,vstep

ALLOCATE(P(vmax-vmin+2))

READ(16,*) numaf
ALLOCATE(airfile(numaf),nAlphaMax(numaf))

DO af=1,numaf
    READ(16,*) airfile(af)
END DO

READ(16,*) imax 
ALLOCATE(rR(imax),cR(imax),tw(imax),ai(imax), &
   Vn(imax),Vt(imax),phi(imax),beta(imax),iAlpha(imax), &
   cli(imax),cdi(imax),Lp(imax),Dp(imax),Psec(imax),phir(i),Cnor(imax),Ctan(imax))

DO i=1,imax
    READ(16,*) rR(i),cR(i),tw(i),ai(i)
END DO

! THIS FILE NEEDS TO BE IN CURRENT FOLDER, FORTRAN WONT READ "/", FIXABLE??? (below)
DO af=1,numaf
    OPEN(af*10,FILE = airfile(af), FORM = 'FORMATTED')
    READ(af*10,*) nAlphaMax(af)
    ALLOCATE(alpha(nAlphaMax(af),numaf),cl(nAlphaMax(af),numaf),cd(nAlphaMax(af),numaf))
    DO nAlpha = 1,nAlphaMax(af)
        READ(af*10,*) alpha(nAlpha,af),cl(nAlpha,af),cd(nAlpha,af) 
    END DO
END DO

omega = rpm*2.0*pi/60.0
wSec = (R-rootcut*R)/imax
nmax = 20         !<=========== max number of iterations (same 'problem' as griddgen)
ALLOCATE(a(nmax),ap(nmax))
!------------------------------------------------------------------------------------
!===================== MAIN LOOP STARTS HERE ========================================
!------------------------------------------------------------------------------------
DO V = vmin,vmax,vstep
    DO i = 1,imax
        !--------------------- Initialize a and ap -------------------------------------------
        lam = omega*rR(i)*R/V
        sig = (B*cR(i)*R)/(2*pi*rR(i)*R)
        beta(i) = tw(i) + pitch
        a(1) = 0.25*(2.0+(pi*lam*sig)-SQRT(4-(4*pi*lam*sig)+(pi*(lam**2)*sig) &
            *((8*pi*beta(i)/180.0)+(pi*sig))))
        ap(1) = 0
        !-------------------- iteratively solve for a and ap ---------------------------------  
        DO n = 1,nmax
            Vn(i) = V*(1-a(n))
            Vt(i) = omega*rR(i)*R*(1+ap(n))
            phi(i) = 180*(ATAN2(Vn(i),Vt(i)))/pi
            beta(i) = tw(i) + pitch
            iAlpha(i) = phi(i) - beta(i)
           ! WRITE(6,*) iAlpha(i)

            !------------ interpolate airfoil file for cl and cd --------------------------------
            af = ai(i)
            DO nAlpha = 1,nAlphaMax(af)-1
                IF (alpha(nAlpha+1,af) > iAlpha(i) .AND. alpha(nAlpha,af) < iAlpha(i)) THEN

                    cli(i) = cl(nAlpha,af)+((cl(nAlpha+1,af)-cl(nAlpha,af))/ &
                        (alpha(nAlpha+1,af)-alpha(nAlpha,af)))*(iAlpha(i)-alpha(nAlpha,af))

                    cdi(i) = cd(nAlpha,af)+((cd(nAlpha+1,af)-cd(nAlpha,af))/ &
                        (alpha(nAlpha+1,af)-alpha(nAlpha,af)))*(iAlpha(i)-alpha(nAlpha,af))
                    EXIT
                END IF
            END DO

            !------------- Redetermine a and ap (with cl and cd) ---------------------------------
            phir(i) = pi*phi(i)/180 
            F = 2.0*ACOS(EXP(-(B*R*(1-rR(i)))/(2*R*rR(i)*SIN(phir(i)))))/pi
            A1 = 4.0*F*((SIN(phir(i)))**2)
            B1 = sig*(cli(i)*COS(phir(i))+cdi(i)*SIN(phir(i)))

            a(n+1) = (1.0+A1/b1)**(-1.0)
            ap(n+1) = (-1.0+(4.0*F *SIN(phir(i))*COS(phir(i)))/ &
            (sig*(cli(i)*SIN(phir(i))-cdi(i)*COS(phir(i)))))**(-1.0)

            IF( (ABS(a(n+1)-a(n))+ABS(ap(n+1)-ap(n))) <= 1.0e-9) THEN 
                EXIT
            END IF
        END DO

        !------------- Find Sectional lift, drag, and power ---------------------------------- 
        Lp(i) = 0.5*rho*(Vn(i)**2 + Vt(i)**2)*cR(i)*R*cli(i)
        Dp(i) = 0.5*rho*(Vn(i)**2 + Vt(i)**2)*cR(i)*R*cdi(i)

        Ctan(i) = cli(i)*SIN(pi*iAlpha(i)/180) - cdi(i)*COS(pi*iAlpha(i)/180)
        Cnor(i) = cli(i)*COS(pi*iAlpha(i)/180) + cdi(i)*SIN(pi*iAlpha(i)/180)

        IF (V == 7) THEN
        WRITE(6,*) rR(i),Ctan(i),Cnor(i)
        END IF 

        Psec(i) = (Lp(i)*SIN(Phir(i))-Dp(i)*COS(Phir(i)))*rR(i)*R*omega*wSec
        P(V) = P(V) + B*Psec(i)

    END DO
    WRITE(6,*) v,P(V)
END DO  

END PROGRAM wT_ae6701