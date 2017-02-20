PROGRAM elliptic_grid_gen
    IMPLICIT NONE
    !
    !----------------------------Variables ---------------------------------------------------
    INTEGER,PARAMETER:: rDef=SELECTED_REAL_KIND(12)
    INTEGER::i,j,n,k,p,imax,jmax,nmax,itr
    REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:,:) :: X,Xsi,Xsisi,Xeta,Xetaeta, &
    Y,Ysi,Ysisi,Yeta,Yetaeta, &
    alpha,beta,gama, &
    bbx,ddx,aax,ccx, &
    bby,ddy,aay,ccy,xResid,yResid

    REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:) :: psi,phi
    REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:) :: RMSres,Sres
    CHARACTER(len=30) :: infile,outfile
    CHARACTER(len=8) :: junk
    LOGICAL :: srctrms
    REAL(KIND=rDef) :: cRMSr
    !----------------------- Open files and set stipulations --------------------------------------

    OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
    READ(16,*) infile 
    READ(16,*) outfile
    READ(16,*) srctrms
    READ(16,*) nmax
    READ(16,*) cRMSr

    OPEN(26,FILE = infile, FORM = 'FORMATTED')

    READ(26,'(A)') junk 
    READ(26,'(A,I3,A,I3)') junk,imax,junk,jmax
    
    ALLOCATE(X(imax,jmax,nmax+1),Xsi(imax,jmax,nmax),Xsisi(imax,jmax,nmax),Xeta(imax,jmax,nmax), &
        Xetaeta(imax,jmax,nmax), &
        Y(imax,jmax,nmax+1),Ysi(imax,jmax,nmax),Ysisi(imax,jmax,nmax),Yeta(imax,jmax,nmax), &
        Yetaeta(imax,jmax,nmax), &
        alpha(imax,jmax,nmax),beta(imax,jmax,nmax),gama(imax,jmax,nmax), &
        bbx(imax,jmax,nmax),ddx(imax,jmax,nmax),aax(imax,jmax,nmax),ccx(imax,jmax,nmax), &
        bby(imax,jmax,nmax),ddy(imax,jmax,nmax),aay(imax,jmax,nmax),ccy(imax,jmax,nmax), &
        xResid(imax,jmax,nmax+1),yResid(imax,jmax,nmax+1))
    ALLOCATE(psi(imax,jmax),phi(imax,jmax)) 
    ALLOCATE(RMSres(nmax),Sres(nmax))

    n = 1                                                            
    DO j=1,jmax                                                      
        DO i=1,imax                                                  
            READ(26,*) X(i,j,n),Y(i,j,n)
        END DO 
    END DO
    CLOSE(26) 

    !-------------------- Find Phi and Psi (constant with n; woohoo!) ------------------------
    n = 1
    DO j=1,jmax,(jmax-1)
        DO i=2,imax-1
            IF (ABS(X(i+1,j,n)-X(i-1,j,n)) > ABS(Y(i+1,j,n)-Y(i-1,j,n))) THEN
                phi(i,j) = -2*(X(i+1,j,n)-2*X(i,j,n)+X(i-1,j,n))/(X(i+1,j,n)-X(i-1,j,n))
            ELSE
                phi(i,j) = 2*(Y(i+1,j,n)-2*Y(i,j,n)+Y(i-1,j,n))/(Y(i+1,j,n)-Y(i-1,j,n))
            END IF
        END DO
    END DO

    DO i=1,imax,(imax-1)
        DO j=2,jmax-1
            IF (ABS(X(i,j+1,n)-X(i,j-1,n)) > ABS(Y(i,j+1,n)-Y(i,j-1,n))) THEN
                psi(i,j) = -2*(X(i,j+1,n)-2*X(i,j,n)+X(i,j-1,n))/(X(i,j+1,n)-X(i,j-1,n))
            ELSE
                psi(i,j) = -2*(Y(i,j+1,n)-2*Y(i,j,n)+Y(i,j-1,n))/(Y(i,j+1,n)-Y(i,j-1,n))
            END IF
        END DO
    END DO

    DO j=2,jmax-1
        DO i=2,imax
            phi(i,j) = phi(i,1)+((j-1.0)/(jmax-1.0))*(phi(i,jmax)-phi(i,1))
            psi(i,j) = psi(1,j)+((i-1.0)/(imax-1.0))*(psi(imax,j)-psi(1,j))
        END DO
    END DO
    

!--------------------------------------------------------------------------------------------------------
!======================= MAIN LOOP STARTS HERE ==========================================================
!--------------------------------------------------------------------------------------------------------
DO n = 1,nmax
    Sres(n) = 0
    !----------------------- set boundaries to be constant -----------------------------------
    DO i=1,imax
        X(i,1,n+1) = X(i,1,1)               
        X(i,jmax,n+1) = X(i,jmax,1)           
        Y(i,1,n+1) = Y(i,1,1)
        Y(i,jmax,n+1) = Y(i,jmax,1)          
    END DO 

    DO j=1,jmax
        X(1,j,n+1) = X(1,j,1)
        X(imax,j,n+1) = X(imax,j,1)
        Y(1,j,n+1) = Y(1,j,1)
        Y(imax,j,n+1) = Y(imax,j,1)
    END DO

    DO j=2,jmax-1
        !-------------------------- get alpha, beta, and gama ---------------------------------------
        DO i=2,imax-1
            Xeta(i,j,n) = (X(i,j+1,n)-X(i,j-1,n))/2
            Yeta(i,j,n) = (Y(i,j+1,n)-Y(i,j-1,n))/2
            alpha(i,j,n) = ((Xeta(i,j,n))**2)+((Yeta(i,j,n))**2)

            Xsi(i,j,n) = (X(i+1,j,n)-X(i-1,j,n))/2
            Ysi(i,j,n) = (Y(i+1,j,n)-Y(i-1,j,n))/2 
            gama(i,j,n) = ((Xsi(i,j,n))**2)+((Ysi(i,j,n))**2)

            beta(i,j,n) = Xsi(i,j,n)*Xeta(i,j,n) + Ysi(i,j,n)*Yeta(i,j,n)

            !-------------------Finding Thomas array values --------------------------------------------
            IF (srctrms .eqv. .FALSE.) THEN   
                bbx(i,j,n) = alpha(i,j,n)
                aax(i,j,n) = bbx(i,j,n) 
                ccx(i,j,n) = beta(i,j,n)*(X(i+1,j+1,n)-X(i+1,j-1,n+1)-X(i-1,j+1,n)+X(i-1,j-1,n+1))/2 &
                -gama(i,j,n)*(X(i,j+1,n)+X(i,j-1,n+1))

                bby(i,j,n) = alpha(i,j,n)
                aay(i,j,n) = bby(i,j,n) 
                ccy(i,j,n) = beta(i,j,n)*(Y(i+1,j+1,n)-Y(i+1,j-1,n+1)-Y(i-1,j+1,n)+Y(i-1,j-1,n+1))/2 &
                -gama(i,j,n)*(Y(i,j+1,n)+Y(i,j-1,n+1))
            ELSE
                bbx(i,j,n) = alpha(i,j,n)*(1 - phi(i,j)/2)

                aax(i,j,n) = alpha(i,j,n)*(1 + phi(i,j)/2)
                ccx(i,j,n) = beta(i,j,n)*(X(i+1,j+1,n)-X(i+1,j-1,n+1)-X(i-1,j+1,n)+X(i-1,j-1,n+1))/2 &
                -gama(i,j,n)*(X(i,j+1,n)+X(i,j-1,n+1))-gama(i,j,n)*psi(i,j)*(X(i,j+1,n)-X(i,j-1,n+1))/2

                bby(i,j,n) = alpha(i,j,n)*(1 - phi(i,j)/2)
                aay(i,j,n) = alpha(i,j,n)*(1 + phi(i,j)/2)
                ccy(i,j,n) = beta(i,j,n)*(Y(i+1,j+1,n)-Y(i+1,j-1,n+1)-Y(i-1,j+1,n)+Y(i-1,j-1,n+1))/2 &
                -gama(i,j,n)*(Y(i,j+1,n)+Y(i,j-1,n+1))-gama(i,j,n)*psi(i,j)*(Y(i,j+1,n)-Y(i,j-1,n+1))/2

            END IF
            ddy(i,j,n) = (-2)*(alpha(i,j,n) + gama(i,j,n))
            ddx(i,j,n) = (-2)*(alpha(i,j,n) + gama(i,j,n))

            !------------------------ Find innitial Residuals -------------------------------------------
            IF (n == 1) THEN
                xResid(i,j,n) = bbx(i,j,n)*X(i-1,j,n)+ddx(i,j,n)*X(i,j,n)+aax(i,j,n)*X(i+1,j,n)-ccx(i,j,n)
                yResid(i,j,n) = bby(i,j,n)*Y(i-1,j,n)+ddy(i,j,n)*Y(i,j,n)+aay(i,j,n)*Y(i+1,j,n)-ccy(i,j,n)
            END IF
        END DO

        !----------------- Setting Thomas array boundaries -----------------------------------------
        bbx(1,j,n) = 0
        bbx(imax,j,n) = 0
        ddx(1,j,n) = 1
        ddx(imax,j,n) = 1        ! X array boundaries
        aax(1,j,n) = 0
        aax(imax,j,n) = 0
        ccx(1,j,n) = x(1,j,n)
        ccx(imax,j,n) = x(imax,j,n)

        bby(1,j,n) = 0
        bby(imax,j,n) = 0
        ddy(1,j,n) = 1
        ddy(imax,j,n) = 1        ! Y array boundaries
        aay(1,j,n) = 0
        aay(imax,j,n) = 0
        ccy(1,j,n) = y(1,j,n)
        ccy(imax,j,n) = y(imax,j,n)

        !----------------------- Thomas algorithm --------------------------------------
        DO i=2,imax-1
            ddx(i,j,n) = ddx(i,j,n) - (bbx(i,j,n)/ddx(i-1,j,n))*aax(i-1,j,n)  ! X forward sweep
            ccx(i,j,n) = ccx(i,j,n) - (bbx(i,j,n)/ddx(i-1,j,n))*ccx(i-1,j,n)

            ddy(i,j,n) = ddy(i,j,n) - (bby(i,j,n)/ddy(i-1,j,n))*aay(i-1,j,n)  ! Y forward sweep
            ccy(i,j,n) = ccy(i,j,n) - (bby(i,j,n)/ddy(i-1,j,n))*ccy(i-1,j,n)
        END DO

        ccx(imax,j,n) = ccx(imax,j,n)/ddx(imax,j,n)  
        ccy(imax,j,n) = ccy(imax,j,n)/ddy(imax,j,n)  

        DO i=imax,2,-1
            ccx(i,j,n) = (ccx(i,j,n) - aax(i,j,n)*ccx(i+1,j,n))/ddx(i,j,n)   ! X forward sweep
            X(i,j,n+1) = ccx(i,j,n)

            ccy(i,j,n) = (ccy(i,j,n) - aay(i,j,n)*ccy(i+1,j,n))/ddy(i,j,n)   ! Y forward sweep
            Y(i,j,n+1) = ccy(i,j,n)
        END DO

        !----------------------- Find Residual and its sum of squares ------------------------------------
        DO i = 2,imax-1
            xResid(i,j,n+1) = bbx(i,j,n)*X(i-1,j,n)+ddx(i,j,n)*X(i,j,n)+aax(i,j,n)*X(i+1,j,n)-ccx(i,j,n)
            yResid(i,j,n+1) = bby(i,j,n)*Y(i-1,j,n)+ddy(i,j,n)*Y(i,j,n)+aay(i,j,n)*Y(i+1,j,n)-ccy(i,j,n)

            Sres(n) = Sres(n) + (xResid(i,j,n+1) - xResid(i,j,n))**2 + (yResid(i,j,n+1) - yResid(i,j,n))**2
        END DO

    END DO !---->  j loop

    RMSres(n) = SQRT(Sres(n)/((imax-2)*(jmax-2)*2)) 

    IF (RMSres(n) <= cRMSr) THEN
        WRITE(6,*) n,RMSres(n)
        WRITE(6,'(A)') ' '
        WRITE(6,'(A,I3,A)') ' Met convergence criteria in ',n,' iterations'
        WRITE(6,*) 'RMS residual history written above.'
        WRITE(6,*) 'Grid wrote to:   ',outfile
        EXIT
    ELSE
        WRITE(6,*) n,RMSres(n)
    END IF

END DO !---->  n loop (main)
!--------------------------------------------------------------------------------------------------
!======================= END OF MAIN LOOP =========================================================
!--------------------------------------------------------------------------------------------------

    !-----------------------------Just in case -----------------------------------------------------
    IF (n >= nmax .AND. RMSres(nmax) > cRMSr) THEN
        WRITE(6,'(A)') 'CONVERGENCE CRITERIA NOT MET! Grid was still written.'
        WRITE(6,'(A)') 'RMS residual history written above'
    END IF

    !------------------------- Writing Results to file for TecPlot360ex--------------------------------
    OPEN(36,FILE = outfile, FORM = 'FORMATTED')
    WRITE(36,'(A)') 'VARIABLES = "X" "Y" '
    WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
    DO j=1,jmax
      DO i=1,imax
        WRITE(36,*) X(i,j,n),Y(i,j,n)
      END DO
    END DO
    CLOSE(36)

END PROGRAM elliptic_grid_gen