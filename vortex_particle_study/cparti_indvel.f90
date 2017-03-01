PROGRAM particle_induced_v
    IMPLICIT NONE

INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(10)
REAL(KIND=rDef),PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
INTEGER :: nfil,j,NumPpv,n,nppv,p
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:)::xm,ym,zm, &
                                          dxarr,dyarr,dzarr, &
                                          rsarr,rarr

REAL(KIND=rDef):: gama,xpoi,ypoi,zpoi, &
                  xm1,xm2,ym1,ym2,zm1,zm2,xmMag,ymMag,zmMag,mMag, &
                  xSeg,ySeg,zSeg,px,py,pz,rxp,ryp,rzp,rpMag, &
                  dvix,dviy,dviz,gcm,vix,viy,viz, &
                  pxAlph,pyAlph,pzAlph,sigj,rho,G,K,cZeta,Zeta, &
                  xpAlphRp,ypAlphRp,zpAlphRp,Csig,hRes,xstep,xpoii

OPEN(16,FILE = 'vector_and_point.inp', FORM = 'FORMATTED')

READ(16,*) gama
READ(16,*) Csig
READ(16,*) hRes
READ(16,*) xpoi,ypoi,zpoi
READ(16,*) NumPpv
READ(16,*) nfil
ALLOCATE(xm(nfil),ym(nfil),zm(nfil), &
         dxarr(nfil),dyarr(nfil),dzarr(nfil), &
         rsarr(nfil),rarr(nfil))

! IF you want to store 'em (also make them allocatable!)
!ALLOCATE(px((nfil-1)*NumPpv),py((nfil-1)*NumPpv),pz((nfil-1)*NumPpv))

DO j = 1,nfil
    READ(16,*) xm(j),ym(j),zm(j)
END DO 

vix = 0.00
viy = 0.00
viz = 0.00
!-------------------------------------------------------------------------
!======================= MAIN LOOP STARTS HERE =========================== 
!-------------------------------------------------------------------------

WRITE(6,*) '------------------------------------------------------------------'
      WRITE(6,*) 'Partical induced velocity at point P:'
      WRITE(6,*) '  xp                       Vi                     '

xstep = .1
xpoii = .5
DO p = 1,21

vix = 0.00
viy = 0.00
viz = 0.00
Xpoi = xpoii + (p-1)*xstep

j = 1  !DO j = 1,nfil-1 
    xm1 = xm(j)
    ym1 = ym(j)
    zm1 = zm(j)    ! Tail and tip vector points 
    xm2 = xm(j+1)
    ym2 = ym(j+1)
    zm2 = zm(j+1)

    xmMag = (xm2-xm1) 
    ymMag = (ym2-ym1)   ! Component magnitudes
    zmMag = (zm2-zm1)
    mMag = SQRT(xmMag**2 + ymMag**2 + zmMag**2)

    ! Defined as "constants" in "n" loop below:
    xSeg = xmMag/(NumPpv*2)
    ySeg = ymMag/(NumPpv*2)
    zSeg = zmMag/(NumPpv*2) 

    pxAlph = gama*xmMag/(NumPpv) !--- Eqv to alpha vector components
    pyAlph = gama*ymMag/(NumPpv) !--- Do these strengths work????
    pzAlph = gama*zmMag/(NumPpv) !            (probably)

    DO n = 1,NumPpv 
       !nppv = n+(j-1)*NumPpv
       px = xm1 + (2*n-1)*xSeg
       py = ym1 + (2*n-1)*ySeg
       pz = zm1 + (2*n-1)*zSeg

       !------ Write out partical point location (if you want) 
       ! WRITE(6,*) nppv,px,py,pz
       rxp = xpoi - px 
       ryp = ypoi - py
       rzp = zpoi - pz

       rpMag = SQRT(rxp**2 + ryp**2 + rzp**2)
        
       ! Find the commonly used Rho (SUPER IMPORTANT)
       sigj = Csig*hRes
       rho = rpMag/sigj

       ! Get K(rho) the kernal used-------
       G = (ERF(rho/root2))/(4*pi*rho) 
       cZeta = 1.0/((2.0*pi)**(3.0/2.0))
       Zeta = cZeta*exp(-((rho)**2.0)/2.0)
       !----------------------------------
       K = (G - Zeta)/(rho**2.0)

       ! Cross product of pAlph and rp (pAlph X rp)
       xpAlphRp  =  pyAlph*rzp - pzAlph*ryp
       ypAlphRp  =  pzAlph*rxp - pxAlph*rzp
       zpAlphRp  =  pxAlph*ryp - pyAlph*rxp

       dvix = xpAlphRp*K/(sigj**3.0)
       dviy = ypAlphRp*K/(sigj**3.0)
       dviz = zpAlphRp*K/(sigj**3.0)

       vix = vix + dvix
       viy = viy + dviy
       viz = viz + dviz

    END DO
!END DO
      WRITE(6,*) xpoi,viz
 END DO 


END PROGRAM particle_induced_v