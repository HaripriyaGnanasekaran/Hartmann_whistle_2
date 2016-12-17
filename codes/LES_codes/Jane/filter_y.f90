        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 26-10-2009        Last Modified on 10-05-2010

        !______________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE FILTERING_Y(u_cur)

        USE GRID_DIMENSIONS,     ONLY: NXL,NXR,NYB,NYT,NZB,NZF
        USE FILTERING_PARAMETERS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE  :: u_fil
        REAL (KIND=8), DIMENSION(NYB:NYT-4,NXL:NXR,NZB:NZF), SAVE   :: FY

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
        do i = NXL,NXR

        u_fil(i,NYB,k)  = a1*u_cur(i,NYB,k) + a2*u_cur(i,NYB+1,k) + a3*u_cur(i,NYB+2,k)
        u_fil(i,NYB,k)  = u_fil(i,NYB,k) + a4*u_cur(i,NYB+3,k) + a5*u_cur(i,NYB+4,k)

        u_fil(i,NYB+1,k)  = b1*u_cur(i,NYB,k) + b2*u_cur(i,NYB+1,k) + b3*u_cur(i,NYB+2,k)
        u_fil(i,NYB+1,k)  = u_fil(i,NYB+1,k) + b4*u_cur(i,NYB+3,k) + b5*u_cur(i,NYB+4,k)

        u_fil(i,NYT,k)   = a1*u_cur(i,NYT,k) + a2*u_cur(i,NYT-1,k)
        u_fil(i,NYT,k)   = u_fil(i,NYT,k) + a3*u_cur(i,NYT-2,k) + a4*u_cur(i,NYT-3,k) 
        u_fil(i,NYT,k)   = u_fil(i,NYT,k) + a5*u_cur(i,NYT-4,k)

        u_fil(i,NYT-1,k) = b1*u_cur(i,NYT,k) + b2*u_cur(i,NYT-1,k)
        u_fil(i,NYT-1,k) = u_fil(i,NYT-1,k) + b3*u_cur(i,NYT-2,k) + b4*u_cur(i,NYT-3,k)
        u_fil(i,NYT-1,k) = u_fil(i,NYT-1,k) + b5*u_cur(i,NYT-4,k)
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL SHARED(u_cur,FY) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)         
        do k = NZB,NZF
        do i = NXL,NXR
        do j = NYB+2,NYT-2
        FY(j-2,i,k) = ay(i)*u_cur(i,j,k) + by(i)*(u_cur(i,j-1,k) + u_cur(i,j+1,k)) 
        FY(j-2,i,k) = FY(j-2,i,k) + cy(i)*(u_cur(i,j+2,k) + u_cur(i,j-2,k))
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_fil,FY) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)         
        do k = NZB,NZF
        do i = NXL,NXR
        FY(NYB,i,k)   = FY(NYB,i,k) - alphay(i)*u_fil(i,NYB+1,k)
        FY(NYT-4,i,k) = FY(NYT-4,i,k) - alphay(i)*u_fil(i,NYT-1,k)
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(FY) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)         
        do k = NZB,NZF
        do i = NXL,NXR
        do j = NYB,NYT-5
        FY(j+1,i,k) = FY(j+1,i,k) - inv_y(i,j)*FY(j,i,k)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_fil,FY) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)         
        do k = NZB,NZF
        do i = NXL,NXR        
        u_fil(i,NYT-2,k) = FY(NYT-4,i,k)*A_coef_y(i,NYT-4,2)
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_fil,FY) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)         
        do k = NZB,NZF
        do i = NXL,NXR     
        do j = NYT-3,NYB+2,-1
        u_fil(i,j,k) = FY(j-2,i,k) - A_coef_y(i,j-2,3)*u_fil(i,j+1,k)
        u_fil(i,j,k) = u_fil(i,j,k)*A_coef_y(i,j-2,2)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL



        CALL COPY(u_fil,u_cur)

        END SUBROUTINE
        !______________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE INITIALIZE_FILTERING_COEFF_Y
        
        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dy1,dz1,jet_r,jet_d
        USE FILTERING_PARAMETERS

        IMPLICIT NONE

        INTEGER :: i,j,k,l
        REAL (KIND = 8) :: x,y,z,r
        
        do i = NXL,NXR

        CALL X_COORD(x,i)
        !alphay(i) = 0.4865d0 - 0.0115*tanh((x-0.04d0)/(4.0d0*jet_r/3.0d0))
        alphay(i) = 0.498d0

        ay(i) = (5.0d0 + 6.0d0*alphay(i))*0.125d0          !/8.0d0
        by(i) = (1.0d0 + 2.0d0*alphay(i))*0.25d0           !/2.0d0*0.5d0
        cy(i) = 0.0d0 - (1.0d0 - 2.0d0*alphay(i))*0.0625d0 !/8.0d0*0.5d0

        do j = NYB,NYT-4
        A_coef_y(i,j,1) = alphay(i)
        A_coef_y(i,j,2) = 1.0d0
        A_coef_y(i,j,3) = alphay(i)
        end do

        A_coef_y(i,NYB,1) = 0.0d0
        A_coef_y(i,NYT-4,3) = 0.0d0

        do j = NYB+1,NYT-4

        inv_y(i,j-1) = A_coef_y(i,j,1)/A_coef_y(i,j-1,2)

        A_coef_y(i,j,1) = 0.0d0
        A_coef_y(i,j,2) = A_coef_y(i,j,2) - A_coef_y(i,j-1,3)*inv_y(i,j-1)
        end do

        do j = NYB,NYT-4
        A_coef_y(i,j,2) = 1.0d0/A_coef_y(i,j,2)
        end do

        
        end do

        END SUBROUTINE
        !______________________________________________________________________
