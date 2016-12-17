        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 26-10-2009        Last Modified on 10-05-2010

        !______________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE FILTERING_Z(u_cur)

        USE GRID_DIMENSIONS,     ONLY: NXL,NXR,NYB,NYT,NZB,NZF
        USE FILTERING_PARAMETERS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE  :: u_fil
        REAL (KIND=8), DIMENSION(NZB:NZF-4,NXL:NXR,NYB:NYT), SAVE   :: FZ

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



!$OMP PARALLEL SHARED(u_cur,u_fil) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)         
        do j = NYB,NYT
        do i = NXL,NXR

        u_fil(i,j,NZB)  = a1*u_cur(i,j,NZB) + a2*u_cur(i,j,NZB+1) + a3*u_cur(i,j,NZB+2)
        u_fil(i,j,NZB)  = u_fil(i,j,NZB) + a4*u_cur(i,j,NZB+3) + a5*u_cur(i,j,NZB+4)

        u_fil(i,j,NZB+1)  = b1*u_cur(i,j,NZB) + b2*u_cur(i,j,NZB+1) + b3*u_cur(i,j,NZB+2)
        u_fil(i,j,NZB+1)  = u_fil(i,j,NZB+1) + b4*u_cur(i,j,NZB+3) + b5*u_cur(i,j,NZB+4)

        u_fil(i,j,NZF)   = a1*u_cur(i,j,NZF) + a2*u_cur(i,j,NZF-1)
        u_fil(i,j,NZF)   = u_fil(i,j,NZF) + a3*u_cur(i,j,NZF-2) + a4*u_cur(i,j,NZF-3) 
        u_fil(i,j,NZF)   = u_fil(i,j,NZF) + a5*u_cur(i,j,NZF-4)

        u_fil(i,j,NZF-1) = b1*u_cur(i,j,NZF) + b2*u_cur(i,j,NZF-1)
        u_fil(i,j,NZF-1) = u_fil(i,j,NZF-1) + b3*u_cur(i,j,NZF-2) + b4*u_cur(i,j,NZF-3)
        u_fil(i,j,NZF-1) = u_fil(i,j,NZF-1) + b5*u_cur(i,j,NZF-4)
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


!$OMP PARALLEL SHARED(u_cur,FZ) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)         
        do j = NYB,NYT
        do i = NXL,NXR
        do k = NZB+2,NZF-2
        FZ(k-2,i,j) = az(i)*u_cur(i,j,k) + bz(i)*(u_cur(i,j,k-1) + u_cur(i,j,k+1)) 
        FZ(k-2,i,j) = FZ(k-2,i,j) + cz(i)*(u_cur(i,j,k+2) + u_cur(i,j,k-2))
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_fil,FZ) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)         
        do j = NYB,NYT
        do i = NXL,NXR
        FZ(NZB,i,j)   = FZ(NZB,i,j) - alphaz(i)*u_fil(i,j,NZB+1)
        FZ(NZF-4,i,j) = FZ(NZF-4,i,j) - alphaz(i)*u_fil(i,j,NZF-1)
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(FZ) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)         
        do j = NYB,NYT
        do i = NXL,NXR
        do k = NZB,NZF-5
        FZ(k+1,i,j) = FZ(k+1,i,j) - inv_z(i,k)*FZ(k,i,j)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_fil,FZ) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)         
        do j = NYB,NYT
        do i = NXL,NXR        
        u_fil(i,j,NZF-2) = FZ(NZF-4,i,j)*A_coef_z(i,NZF-4,2)
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_fil,FZ) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)         
        do j = NYB,NYT
        do i = NXL,NXR     
        do k = NZF-3,NZB+2,-1
        u_fil(i,j,k) = FZ(k-2,i,j) - A_coef_z(i,k-2,3)*u_fil(i,j,k+1)
        u_fil(i,j,k) = u_fil(i,j,k)*A_coef_z(i,k-2,2)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL



        CALL COPY(u_fil,u_cur)

        END SUBROUTINE
        !______________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE INITIALIZE_FILTERING_COEFF_Z
        
        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NZB,NZF
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dy1,dz1,jet_r,jet_d
        USE FILTERING_PARAMETERS

        IMPLICIT NONE

        INTEGER :: i,j,k,l
        REAL (KIND = 8) :: x,y,z,r

        do i = NXL,NXR

        CALL X_COORD(x,i)
        !alphaz(i) = 0.4865d0 - 0.0115*tanh((x-0.04d0)/(4.0d0*jet_r/3.0d0))        
        alphaz(i) = 0.498d0

        az(i) = (5.0d0 + 6.0d0*alphaz(i))*0.125d0          !/8.0d0
        bz(i) = (1.0d0 + 2.0d0*alphaz(i))*0.25d0           !/2.0d0*0.5d0
        cz(i) = 0.0d0 - (1.0d0 - 2.0d0*alphaz(i))*0.0625d0 !/8.0d0*0.5d0
        
        do k = NZB,NZF-4
        A_coef_z(i,k,1) = alphaz(i)
        A_coef_z(i,k,2) = 1.0d0
        A_coef_z(i,k,3) = alphaz(i)
        end do

        A_coef_z(i,NZB,1) = 0.0d0
        A_coef_z(i,NZF-4,3) = 0.0d0

        do k = NZB+1,NZF-4

        inv_z(i,k-1) = A_coef_z(i,k,1)/A_coef_z(i,k-1,2)

        A_coef_z(i,k,1) = 0.0d0
        A_coef_z(i,k,2) = A_coef_z(i,k,2) - A_coef_z(i,k-1,3)*inv_z(i,k-1)
        end do

        do k = NZB,NZF-4
        A_coef_z(i,k,2) = 1.0d0/A_coef_z(i,k,2)
        end do

    
        end do


        END SUBROUTINE
        !______________________________________________________________________
