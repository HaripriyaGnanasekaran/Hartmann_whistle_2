        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 26-10-2009        Last Modified on 10-05-2010

        !______________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE FILTERING_X(u_cur)

        USE GRID_DIMENSIONS,     ONLY: NXL,NXR,NYB,NYT,NZB,NZF
        USE FILTERING_PARAMETERS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE  :: u_fil
        REAL (KIND=8), DIMENSION(NXL:NXR-4,NYB:NYT,NZB:NZF), SAVE   :: FX

        REAL (KIND=8), PARAMETER :: a1 = 15.0d0*0.0625d0
        REAL (KIND=8), PARAMETER :: a2 = 4.0d0*0.0625d0
        REAL (KIND=8), PARAMETER :: a3 = (0.0d0 - 6.0d0)*0.0625d0
        REAL (KIND=8), PARAMETER :: a4 = 4.0d0*0.0625d0
        REAL (KIND=8), PARAMETER :: a5 = (0.0d0 - 1.0d0)*0.0625d0

        REAL (KIND=8), PARAMETER :: b2 = 12.0d0*0.0625d0
        REAL (KIND=8), PARAMETER :: b1 = 1.0d0*0.0625d0
        REAL (KIND=8), PARAMETER :: b3 = 6.0d0*0.0625d0
        REAL (KIND=8), PARAMETER :: b4 = (0.0d0 - 4.0d0)*0.0625d0
        REAL (KIND=8), PARAMETER :: b5 = 1.0d0*0.0625d0


        INTEGER :: i,j,k



!$OMP PARALLEL SHARED(u_cur,u_fil) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)        
        do k = NZB,NZF
        do j = NYB,NYT

        u_fil(NXL,j,k) = a1*u_cur(NXL,j,k) + a2*u_cur(NXL+1,j,k) + a3*u_cur(NXL+2,j,k)
        u_fil(NXL,j,k) = u_fil(NXL,j,k) + a4*u_cur(NXL+3,j,k) + a5*u_cur(NXL+4,j,k)

        u_fil(NXL+1,j,k) = b1*u_cur(NXL,j,k) + b2*u_cur(NXL+1,j,k) + b3*u_cur(NXL+2,j,k)
        u_fil(NXL+1,j,k) = u_fil(NXL+1,j,k) + b4*u_cur(NXL+3,j,k) + b5*u_cur(NXL+4,j,k)

        u_fil(NXR,j,k) = a1*u_cur(NXR,j,k) + a2*u_cur(NXR-1,j,k)
        u_fil(NXR,j,k) = u_fil(NXR,j,k) + a3*u_cur(NXR-2,j,k) + a4*u_cur(NXR-3,j,k) 
        u_fil(NXR,j,k) = u_fil(NXR,j,k) + a5*u_cur(NXR-4,j,k)

        u_fil(NXR-1,j,k) = b1*u_cur(NXR,j,k) + b2*u_cur(NXR-1,j,k)
        u_fil(NXR-1,j,k) = u_fil(NXR-1,j,k) + b3*u_cur(NXR-2,j,k) + b4*u_cur(NXR-3,j,k)
        u_fil(NXR-1,j,k) = u_fil(NXR-1,j,k) + b5*u_cur(NXR-4,j,k)
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_cur,FX) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)        
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL+2,NXR-2
        FX(i-2,j,k) = ax*u_cur(i,j,k) + bx*(u_cur(i-1,j,k) + u_cur(i+1,j,k))
        FX(i-2,j,k) = FX(i-2,j,k) + cx*(u_cur(i+2,j,k) + u_cur(i-2,j,k))
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_fil,FX) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)        
        do k = NZB,NZF
        do j = NYB,NYT
        FX(NXL,j,k)   = FX(NXL,j,k) - alphax*u_fil(NXL+1,j,k)
        FX(NXR-4,j,k) = FX(NXR-4,j,k) - alphax*u_fil(NXR-1,j,k)
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(FX) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)        
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR-5
        FX(i+1,j,k) = FX(i+1,j,k) - inv_x(i)*FX(i,j,k)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_fil,FX) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)        
        do k = NZB,NZF
        do j = NYB,NYT        
        u_fil(NXR-2,j,k) = FX(NXR-4,j,k)*A_coef_x(NXR-4,2)
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_fil,FX) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)        
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXR-3,NXL+2,-1
        u_fil(i,j,k) = FX(i-2,j,k) - A_coef_x(i-2,3)*u_fil(i+1,j,k)
        u_fil(i,j,k) = u_fil(i,j,k)*A_coef_x(i-2,2)
        end do
        end do 
        end do
!$OMP END DO
!$OMP END PARALLEL


        CALL COPY(u_fil,u_cur)

        END SUBROUTINE
        !______________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE INITIALIZE_FILTERING_COEFF_X
        
        USE GRID_DIMENSIONS, ONLY : NXL,NXR
        USE FILTERING_PARAMETERS

        IMPLICIT NONE

        INTEGER :: i,j
        
        
        do i = NXL,NXR-4
        A_coef_x(i,1) = alphax
        A_coef_x(i,2) = 1.0d0
        A_coef_x(i,3) = alphax
        end do

        A_coef_x(NXL,1) = 0.0d0
        A_coef_x(NXR-4,3) = 0.0d0


        do i = NXL+1,NXR-4

        inv_x(i-1) = A_coef_x(i,1)/A_coef_x(i-1,2)

        A_coef_x(i,1) = 0.0d0
        A_coef_x(i,2) = A_coef_x(i,2) - A_coef_x(i-1,3)*inv_x(i-1)
        end do


        do i = NXL,NXR-4
        A_coef_x(i,2) = 1.0d0/A_coef_x(i,2)
        end do

        END SUBROUTINE
        !______________________________________________________________________
