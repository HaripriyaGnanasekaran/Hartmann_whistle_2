
        !______________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 04-05-2010        Last Modified on 04-05-2010

        !______________________________________________________________________





        !______________________________________________________________________
        SUBROUTINE PRODUCT_OF_TWO_TERMS(one_cur,two_cur,res_cur)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        IMPLICIT NONE
        
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)  :: one_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)  :: two_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(OUT) :: res_cur
       
        INTEGER :: i,j,k
        
!$OMP PARALLEL SHARED(one_cur,two_cur,res_cur) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT        
        do i = NXL,NXR
        res_cur(i,j,k) =  one_cur(i,j,k) * two_cur(i,j,k)
        end do
        end do     
        end do
!$OMP END DO
!$OMP END PARALLEL

        END SUBROUTINE
        !______________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE MULTIPLY(one_cur,res_cur)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        IMPLICIT NONE
        
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)  :: one_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(OUT) :: res_cur
       
        INTEGER :: i,j,k
        
!$OMP PARALLEL SHARED(one_cur,res_cur) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        res_cur(i,j,k) =  one_cur(i,j,k) * res_cur(i,j,k)
        end do
        end do
        end do     
!$OMP END DO
!$OMP END PARALLEL

        END SUBROUTINE
        !______________________________________________________________________









        !______________________________________________________________________
        SUBROUTINE FACTOR_OF_TWO_TERMS(one_cur,two_cur,res_cur)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        IMPLICIT NONE
        
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)  :: one_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)  :: two_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(OUT) :: res_cur
       
        INTEGER :: i,j,k
        
!$OMP PARALLEL SHARED(one_cur,two_cur,res_cur) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        res_cur(i,j,k) =  one_cur(i,j,k)/two_cur(i,j,k)
        end do
        end do
        end do     
!$OMP END DO
!$OMP END PARALLEL

        END SUBROUTINE
        !______________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE DIVIDE(one_cur,res_cur)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        IMPLICIT NONE
        
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)  :: one_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(OUT) :: res_cur
       
        INTEGER :: i,j,k
        
!$OMP PARALLEL SHARED(one_cur,res_cur) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        res_cur(i,j,k) =  res_cur(i,j,k)/one_cur(i,j,k)
        end do
        end do
        end do     
!$OMP END DO
!$OMP END PARALLEL

        END SUBROUTINE
        !______________________________________________________________________





        !______________________________________________________________________
        SUBROUTINE SUM_OF_TWO_TERMS(one_cur,two_cur,res_cur)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        IMPLICIT NONE
        
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)  :: one_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)  :: two_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(OUT) :: res_cur
       
        INTEGER :: i,j,k
        
!$OMP PARALLEL SHARED(one_cur,two_cur,res_cur) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        res_cur(i,j,k) =  one_cur(i,j,k) + two_cur(i,j,k)
        end do
        end do
        end do     
!$OMP END DO
!$OMP END PARALLEL

        END SUBROUTINE
        !______________________________________________________________________





        !______________________________________________________________________
        SUBROUTINE ADD(one_cur,res_cur)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        IMPLICIT NONE
        
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)  :: one_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(OUT) :: res_cur
       
        INTEGER :: i,j,k
        
!$OMP PARALLEL SHARED(one_cur,res_cur) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        res_cur(i,j,k) =  one_cur(i,j,k) + res_cur(i,j,k)
        end do
        end do
        end do     
!$OMP END DO
!$OMP END PARALLEL

        END SUBROUTINE
        !______________________________________________________________________






        !______________________________________________________________________
        SUBROUTINE X_COORD(x,i)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS, ONLY : NX_st,NXO,NYSB,NYST,NZSB,NZSF
        USE GRID_SIZE_TIME_STEP, ONLY:dx1,dy1,dz1,axo,ayo,azo,axs
        

        IMPLICIT NONE

        REAL (KIND=8),  INTENT(OUT) :: x
        INTEGER,        INTENT(IN)  :: i
        REAL (KIND=8) :: te
        if (i <= NX_st) then
        x = i * dx1
        elseif(i .gt. NX_st .and. i .le. NXO) THEN
        x = NX_st * dx1 + dx1*axs*((axs**(i-NX_st) - 1.0d0)/(axs-1.0d0))
        else
        te = NX_st * dx1 + dx1*axs*((axs**(NXO-NX_st) -1.0d0)/(axs-1.0d0))
        x = te + dx1*axo*((axo**(i-NXO) - 1.0d0)/(axo-1.0d0))
        end if

        END SUBROUTINE
        !______________________________________________________________________









        !______________________________________________________________________
        SUBROUTINE Y_COORD(y,j)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS, ONLY : NXO,NYSB,NYST,NZSB,NZSF
        USE GRID_SIZE_TIME_STEP, ONLY:dx1,dy1,dz1,axo,ayo,azo
        

        IMPLICIT NONE

        REAL (KIND=8),  INTENT(OUT) :: y
        INTEGER,        INTENT(IN)  :: j

        if ((j >= NYSB) .AND. (j <= NYST)) then
        y = j * dy1
        else if ( j > NYST) then
        y = NYST * dy1 + dy1*ayo*((ayo**(j - NYST) - 1.0d0)/(ayo-1.0d0))
        else if (j < NYSB) then
        y = NYSB * dy1 - dy1*ayo*((ayo**(NYSB - j) - 1.0d0)/(ayo-1.0d0))
        end if

        END SUBROUTINE
        !______________________________________________________________________












        !______________________________________________________________________
        SUBROUTINE Z_COORD(z,k)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS, ONLY : NXO,NYSB,NYST,NZSB,NZSF
        USE GRID_SIZE_TIME_STEP, ONLY:dx1,dy1,dz1,axo,ayo,azo
        

        IMPLICIT NONE

        REAL (KIND=8),  INTENT(OUT) :: z
        INTEGER,        INTENT(IN)  :: k

        if ((k >= NZSB) .AND. (k <= NZSF)) then
        z = k * dz1
        else if ( k > NZSF) then
        z = NZSF * dz1 + dz1*azo*((azo**(k - NZSF) - 1.0d0)/(azo-1.0d0))
        else if (k < NZSB) then
        z = NZSB * dz1 - dz1*azo*((azo**(NZSB - k) - 1.0d0)/(azo-1.0d0))
        end if

        END SUBROUTINE
        !______________________________________________________________________








        !______________________________________________________________________
        SUBROUTINE SUTHERLAND(T,mu)

        IMPLICIT NONE

        REAL (KIND=8), INTENT(IN)  :: T
        REAL (KIND=8), INTENT(OUT) :: mu

        REAL (KIND=8) :: c1 = 1.458*10.0d0**(-6)
        REAL (KIND=8) :: c2 = 110.4d0

        mu = c1*T**1.5/(T+c2)

        END SUBROUTINE
        !______________________________________________________________________













