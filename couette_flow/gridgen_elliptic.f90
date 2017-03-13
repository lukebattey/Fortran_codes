PROGRAM elliptic_grid_gen
  IMPLICIT NONE
!----------------------------Variables ----------------------------------------
INTEGER,PARAMETER:: rDef=SELECTED_REAL_KIND(12)
INTEGER::i,j,n,k,p,imax,jmax,nmax,itr
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:) :: X,Y,Ja, &
alpha,beta,gama,bbx,ddx,aax,ccx,bby,ddy,aay,ccy,XNext,YNext,psi,phi
CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk
LOGICAL :: srctrms
REAL(KIND=rDef) :: cRMS,eRMS,XeSumS,YeSumS,Xsi,Xeta,Ysi,Yeta,&
XeRMS,YeRMS,Xes,Yes

!----------------------- Open files and set stipulations ----------------------
OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
READ(16,*) srctrms
READ(16,*) nmax
READ(16,*) cRMS

OPEN(26,FILE = infile, FORM = 'FORMATTED')
READ(26,'(A)') junk 
READ(26,'(A,I3,A,I3)') junk,imax,junk,jmax

ALLOCATE(X(imax,jmax),Y(imax,jmax),Ja(imax,jmax), &
  alpha(imax,jmax),beta(imax,jmax), &
  gama(imax,jmax),bbx(imax,jmax),ddx(imax,jmax),aax(imax,jmax),ccx(imax,jmax),&
  bby(imax,jmax),ddy(imax,jmax),aay(imax,jmax),ccy(imax,jmax), &
  XNext(imax,jmax),YNext(imax,jmax), &
  psi(imax,jmax),phi(imax,jmax)) 

DO j=1,jmax                                                      
  DO i=1,imax                                                  
    READ(26,*) X(i,j),Y(i,j)
  END DO 
END DO
CLOSE(26) 

!-------------- Find Phi and Psi (constant with n; woohoo!) ------------------
n = 1
DO j=1,jmax,(jmax-1)
  DO i=2,imax-1
    IF (ABS(X(i+1,j)-X(i-1,j)) > ABS(Y(i+1,j)-Y(i-1,j))) THEN
      phi(i,j) = -2*(X(i+1,j)-2*X(i,j)+X(i-1,j))/(X(i+1,j)-X(i-1,j))
    ELSE
      phi(i,j) = 2*(Y(i+1,j)-2*Y(i,j)+Y(i-1,j))/(Y(i+1,j)-Y(i-1,j))
    END IF
  END DO
END DO

DO i=1,imax,(imax-1)
  DO j=2,jmax-1
    IF (ABS(X(i,j+1)-X(i,j-1)) > ABS(Y(i,j+1)-Y(i,j-1))) THEN
      psi(i,j) = -2*(X(i,j+1)-2*X(i,j)+X(i,j-1))/(X(i,j+1)-X(i,j-1))
    ELSE
      psi(i,j) = -2*(Y(i,j+1)-2*Y(i,j)+Y(i,j-1))/(Y(i,j+1)-Y(i,j-1))
    END IF
  END DO
END DO

DO j=2,jmax-1
  DO i=2,imax
    phi(i,j) = phi(i,1)+((j-1.0)/(jmax-1.0))*(phi(i,jmax)-phi(i,1))
    psi(i,j) = psi(1,j)+((i-1.0)/(imax-1.0))*(psi(imax,j)-psi(1,j))
  END DO
END DO

!----------------------- set boundaries to be constant ----------------------------------
DO i=1,imax
  XNext(i,1) = X(i,1)               
  XNext(i,jmax) = X(i,jmax)           
  YNext(i,1) = Y(i,1)
  YNext(i,jmax) = Y(i,jmax)          
END DO 

DO j=1,jmax
  XNext(1,j) = X(1,j)
  XNext(imax,j) = X(imax,j)
  YNext(1,j) = Y(1,j)
  YNext(imax,j) = Y(imax,j)
END DO
!----------------------------------------------------------------------------------------
!========================= MAIN LOOP STARTS HERE ========================================
!----------------------------------------------------------------------------------------
DO n = 1,nmax
  XeSumS = 0.0
  YeSumS = 0.0
  DO j=2,jmax-1
    !-------------------------- get alpha, beta, and gama ---------------------------------
    DO i=2,imax-1
      Xeta = (X(i,j+1)-X(i,j-1))/2
      Yeta = (Y(i,j+1)-Y(i,j-1))/2
      alpha(i,j) = ((Xeta)**2)+((Yeta)**2)

      Xsi = (X(i+1,j)-X(i-1,j))/2
      Ysi = (Y(i+1,j)-Y(i-1,j))/2 
      gama(i,j) = ((Xsi)**2)+((Ysi)**2)

      beta(i,j) = Xsi*Xeta + Ysi*Yeta

      !-------------------Finding Thomas array values ---------------------------------
      IF (srctrms .eqv. .FALSE.) THEN   
        bbx(i,j) = alpha(i,j)
        aax(i,j) = bbx(i,j) 
        ccx(i,j) = beta(i,j)*(X(i+1,j+1)-XNext(i+1,j-1)-X(i-1,j+1)+XNext(i-1,j-1))/2 &
        -gama(i,j)*(X(i,j+1)+XNext(i,j-1))

        bby(i,j) = alpha(i,j)
        aay(i,j) = bby(i,j) 
        ccy(i,j) = beta(i,j)*(Y(i+1,j+1)-YNext(i+1,j-1)-Y(i-1,j+1)+YNext(i-1,j-1))/2 &
        -gama(i,j)*(Y(i,j+1)+YNext(i,j-1))
      ELSE
        bbx(i,j) = alpha(i,j)*(1 - phi(i,j)/2)

        aax(i,j) = alpha(i,j)*(1 + phi(i,j)/2)
        ccx(i,j) = beta(i,j)*(X(i+1,j+1)-XNext(i+1,j-1)-X(i-1,j+1)+XNext(i-1,j-1))/2 &
        -gama(i,j)*(X(i,j+1)+XNext(i,j-1))-gama(i,j)*psi(i,j)*(X(i,j+1)-XNext(i,j-1))/2

        bby(i,j) = alpha(i,j)*(1 - phi(i,j)/2)
        aay(i,j) = alpha(i,j)*(1 + phi(i,j)/2)
        ccy(i,j) = beta(i,j)*(Y(i+1,j+1)-YNext(i+1,j-1)-Y(i-1,j+1)+YNext(i-1,j-1))/2 &
        -gama(i,j)*(Y(i,j+1)+YNext(i,j-1))-gama(i,j)*psi(i,j)*(Y(i,j+1)-YNext(i,j-1))/2

      END IF
      ddy(i,j) = (-2)*(alpha(i,j) + gama(i,j))
      ddx(i,j) = (-2)*(alpha(i,j) + gama(i,j))

    END DO

    !----------------- Setting Thomas array boundaries --------------------------------
    bbx(1,j) = 0
    bbx(imax,j) = 0
    ddx(1,j) = 1
    ddx(imax,j) = 1        ! X array boundaries
    aax(1,j) = 0
    aax(imax,j) = 0
    ccx(1,j) = x(1,j)
    ccx(imax,j) = x(imax,j)

    bby(1,j) = 0
    bby(imax,j) = 0
    ddy(1,j) = 1
    ddy(imax,j) = 1        ! Y array boundaries
    aay(1,j) = 0
    aay(imax,j) = 0
    ccy(1,j) = y(1,j)
    ccy(imax,j) = y(imax,j)

    !----------------------- Thomas algorithm --------------------------------------
    DO i=2,imax-1
      ddx(i,j) = ddx(i,j) - (bbx(i,j)/ddx(i-1,j))*aax(i-1,j)  ! X forward sweep
      ccx(i,j) = ccx(i,j) - (bbx(i,j)/ddx(i-1,j))*ccx(i-1,j)

      ddy(i,j) = ddy(i,j) - (bby(i,j)/ddy(i-1,j))*aay(i-1,j)  ! Y forward sweep
      ccy(i,j) = ccy(i,j) - (bby(i,j)/ddy(i-1,j))*ccy(i-1,j)
    END DO

    ccx(imax,j) = ccx(imax,j)/ddx(imax,j)  
    ccy(imax,j) = ccy(imax,j)/ddy(imax,j)  

    DO i=imax,2,-1
      ccx(i,j) = (ccx(i,j) - aax(i,j)*ccx(i+1,j))/ddx(i,j)   ! X forward sweep
      XNext(i,j) = ccx(i,j)

      ccy(i,j) = (ccy(i,j) - aay(i,j)*ccy(i+1,j))/ddy(i,j)   ! Y forward sweep
      YNext(i,j) = ccy(i,j)
    END DO

    !----------------------- Find RMS Error and its sum of squares --------------------
    DO i = 2,imax-1
      Xes = (XNext(i,j) - X(i,j))**2
      Yes = (YNext(i,j) - Y(i,j))**2
      XeSumS = XeSumS + Xes      
      YeSumS = YeSumS + Yes      
    END DO

  END DO !---->  j loop

  eRMS = SQRT((XeSumS+YeSumS)/((imax-2)*(jmax-2)*2))

  XeRMS = SQRT((XeSumS)/((imax-2)*(jmax-2)))
  YeRMS = SQRT((YeSumS)/((imax-2)*(jmax-2)))

  IF (eRMS <= cRMS ) THEN
    WRITE(6,*) n,XeRMS,YeRMS,eRMS
    WRITE(6,'(A)') ' '
    WRITE(6,'(A,I3,A)') 'Met convergence criteria in ',n,' iterations'
    WRITE(6,'(A)') 'RMS error history above: iteration,eRMS,XeRMS,YeRMS'
    WRITE(6,*) 'Grid wrote to:  ',outfile
    EXIT
  ELSE
    WRITE(6,*) n,XeRMS,YeRMS,eRMS
  END IF

  X = XNext
  Y = YNext

END DO !---->  n loop (main)
  !--------------------------------------------------------------------------------------
  !======================= END OF MAIN LOOP =============================================
  !--------------------------------------------------------------------------------------
IF (n >= nmax) THEN  !.AND. eRMS > cRMS
  WRITE(6,'(A)') 'CONVERGENCE CRITERIA NOT MET! Grid wrote to:  ',outfile
  WRITE(6,'(A)') 'RMS error history written above'
END IF

DO j=1,jmax  
  DO i=1,imax
    Xeta = (X(i,j+1)-X(i,j-1))/2
    Yeta = (Y(i,j+1)-Y(i,j-1))/2
    Xsi = (X(i+1,j)-X(i-1,j))/2
    Ysi = (Y(i+1,j)-Y(i-1,j))/2
    IF (j == jmax) THEN
      Xeta = (3*X(i,j)-4*X(i,j-1)+X(i,j-2))/2
      Yeta = (3*Y(i,j)-4*Y(i,j-1)+Y(i,j-2))/2
    END IF
    IF (j == 1) THEN
      Xeta = -1*(3*X(i,j)-4*X(i,j+1)+X(i,j+2))/2
      Yeta = -1*(3*Y(i,j)-4*Y(i,j+1)+Y(i,j+2))/2
    END IF
    IF (i == imax) THEN
      Xsi = (3*X(i,j)-4*X(i-1,j)+X(i-2,j))/2
      Ysi = (3*Y(i,j)-4*Y(i-1,j)+Y(i-2,j))/2
    END IF
    IF (i == 1) THEN
      Xsi = -1*(3*X(i,j)-4*X(i+1,j)+X(i+2,j))/2
      Ysi = -1*(3*Y(i,j)-4*Y(i+1,j)+Y(i+2,j))/2
    END IF
    Ja(i,j) = Xsi*Yeta - Xeta*Ysi
  END DO
END DO

  !--------------------- Writing Results to file for TecPlot360ex------------------------
OPEN(36,FILE = outfile, FORM = 'FORMATTED')
WRITE(36,'(A)') 'VARIABLES = "X" "Y" "J"'
WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
DO j=1,jmax
  DO i=1,imax
    WRITE(36,*) X(i,j),Y(i,j),Ja(i,j)
  END DO
END DO
CLOSE(36)

END PROGRAM elliptic_grid_gen