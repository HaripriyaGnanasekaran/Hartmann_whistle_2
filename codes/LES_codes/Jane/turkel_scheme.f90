        ! Purpose:
        ! This contains the subroutines for calculating the derivatives
        ! using the Hixon method, both forward and backward methods
        !
 
        !______________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 15-10-2008        Last Modified on 12-05-2010

        !______________________________________________________________________


        ! HIXON SCHEMES







        !______________________________________________________________________
        SUBROUTINE TURKEL_BACKDIF_X(u_cur,u_x)
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NXO,NXB,NYSB,NYST,NZSB,NZSF
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dx,axo,dxi,d_eta_x,d_eta_x_inv
        USE HIXON_BACKWARD_PARAMETERS
        USE METRICS,             ONLY : dx_deta
        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),   INTENT(INOUT)    :: u_x
        REAL (KIND=8), DIMENSION(NXL:NXR-1,NYB:NYT,NZB:NZF), SAVE    :: F
        REAL (KIND=8) :: c
        INTEGER :: i,j,k

!$OMP PARALLEL SHARED(u_cur,F) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR-1
        F(i,j,k) = (u_cur(i+1,j,k) - u_cur(i,j,k))*d_eta_x_inv
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL
 
        dx = d_eta_x
     
!$OMP PARALLEL SHARED(u_cur,u_x) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        u_x(NXL,j,k) = 4.0d0*u_cur(NXL+1,j,k) - u_cur(NXL+2,j,k) - 3.0d0*u_cur(NXL,j,k)
        u_x(NXL,j,k) = 0.5d0 * u_x(NXL,j,k)*d_eta_x_inv
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_x,F) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR-1
        u_x(i+1,j,k) = (F(i,j,k) - a*u_x(i,j,k))*adi
	u_x(i+1,j,k) = u_x(i+1,j,k)/dx_deta(i+1)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        !CALL COPY(u_x,u_cur)

        END SUBROUTINE
        !______________________________________________________________________

        !______________________________________________________________________
        SUBROUTINE TURKEL_FORWARD_X(u_cur,u_x)
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NXO,NXB,NYSB,NYST,NZSB,NZSF
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dx,axo,dxi,d_eta_x,d_eta_x_inv
        USE HIXON_BACKWARD_PARAMETERS
        USE METRICS,             ONLY : dx_deta
        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),   INTENT(INOUT)    :: u_x
        REAL (KIND=8), DIMENSION(NXL:NXR-1,NYB:NYT,NZB:NZF), SAVE    :: F
        REAL (KIND=8) :: c
        INTEGER :: i,j,k

!$OMP PARALLEL SHARED(u_cur,F) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR-1
        F(i,j,k) = (u_cur(i+1,j,k) - u_cur(i,j,k))*d_eta_x_inv
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL
 
        dx = d_eta_x
     
!$OMP PARALLEL SHARED(u_cur,u_x) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        u_x(NXR,j,k) = -4.0d0*u_cur(NXR-1,j,k) + u_cur(NXR-2,j,k) + 3.0d0*u_cur(NXR,j,k)
        u_x(NXR,j,k) = 0.5d0 * u_x(NXR,j,k)*d_eta_x_inv
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_x,F) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXR-1,NXL,-1
	u_x(i,j,k) = (F(i,j,k) - a*u_x(i+1,j,k))*adi
      	u_x(i,j,k) = u_x(i,j,k)/dx_deta(i)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        !CALL COPY(u_x,u_cur)

        END SUBROUTINE
        !______________________________________________________________________

        !______________________________________________________________________
        SUBROUTINE TURKEL_BACKDIF_Y(u_cur,u_y)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NYSB,NYST
        USE GRID_SIZE_TIME_STEP, ONLY : dy1,ayo,dy,dyi,d_eta_y,d_eta_y_inv
        USE HIXON_BACKWARD_PARAMETERS
        USE METRICS,             ONLY : dy_deta

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),   INTENT(INOUT) :: u_y
        REAL (KIND=8), DIMENSION(NYB:NYT-1,NXL:NXR,NZB:NZF), SAVE :: F

        REAL (KIND=8) :: c
        INTEGER :: i,j,k

       
!$OMP PARALLEL SHARED(u_cur,F) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        do j = NYB,NYT-1
        F(j,i,k) = (u_cur(i,j+1,k) - u_cur(i,j,k))*d_eta_y_inv
        !dy = dy/ayo
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

      
!$OMP PARALLEL SHARED(u_cur,u_y) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        u_y(i,NYB,k) = 4.0d0*u_cur(i,NYB+1,k) - u_cur(i,NYB+2,k) -3.0d0*u_cur(i,NYB,k)
        u_y(i,NYB,k) = 0.5d0 * u_y(i,NYB,k)*d_eta_y_inv
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(F,u_y) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        do j = NYB,NYT-1
        u_y(i,j+1,k) = (F(j,i,k) - a*u_y(i,j,k))*adi
        u_y(i,j+1,k) = u_y(i,j+1,k)/dy_deta(j+1)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        !CALL COPY(u_y,u_cur)

        END SUBROUTINE
        !______________________________________________________________________

















        !______________________________________________________________________
        SUBROUTINE TURKEL_FORWARD_Y(u_cur,u_y)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NYSB,NYST
        USE GRID_SIZE_TIME_STEP, ONLY : dy1,ayo,dy,dyi,d_eta_y,d_eta_y_inv
        USE HIXON_FORWARD_PARAMETERS
        USE METRICS,             ONLY : dy_deta
        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),   INTENT(INOUT) :: u_y
        REAL (KIND=8), DIMENSION(NYB:NYT-1,NXL:NXR,NZB:NZF), SAVE :: F

        REAL (KIND=8) :: c
        INTEGER :: i,j,k

     
!$OMP PARALLEL SHARED(u_cur,F) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        do j = NYB,NYT-1
        F(j,i,k) = (u_cur(i,j+1,k) - u_cur(i,j,k))*d_eta_y_inv
        !dy = dy/ayo
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

       
!$OMP PARALLEL SHARED(u_cur,u_y) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        u_y(i,NYT,k) = -4.0d0*u_cur(i,NYT-1,k) + u_cur(i,NYT-2,k) + 3.0d0*u_cur(i,NYT,k)
        u_y(i,NYT,k) = 0.5d0 * u_y(i,NYT,k)*d_eta_y_inv
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

        j = NYT-1

!$OMP PARALLEL SHARED(F,u_y) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXl,NXR
        do j = NYT-1,NYB,-1
        u_y(i,j,k) = (F(j,i,k) - a*u_y(i,j+1,k))*adi
        u_y(i,j,k) = u_y(i,j,k)/dy_deta(j)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        !CALL COPY(u_y,u_cur)

        END SUBROUTINE
        !______________________________________________________________________











        !______________________________________________________________________
        SUBROUTINE TURKEL_BACKDIF_Z(u_cur,u_z)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NZSB,NZSF
        USE GRID_SIZE_TIME_STEP, ONLY : dz1,azo,dz,dzi,d_eta_z,d_eta_z_inv
        USE HIXON_BACKWARD_PARAMETERS
        USE METRICS
        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),   INTENT(INOUT) :: u_z
        REAL (KIND=8), DIMENSION(NZB:NZF-1,NXL:NXR,NYB:NYT), SAVE :: F

        REAL (KIND=8) :: c
        INTEGER :: i,j,k


!$OMP PARALLEL SHARED(u_cur,F) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        do k = NZB,NZF-1
        F(k,i,j) = (u_cur(i,j,k+1) - u_cur(i,j,k))*d_eta_z_inv
        !dz = dz/azo
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

  
!$OMP PARALLEL SHARED(u_cur,u_z) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        u_z(i,j,NZB) = 4.0d0*u_cur(i,J,NZB+1) - u_cur(i,j,NZB+2) -3.0d0*u_cur(i,j,NZB)
        u_z(i,j,NZB) = 0.5d0 * u_z(i,j,NZB)*d_eta_z_inv
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(F,u_z) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        do k = NZB,NZF-1
        u_z(i,j,k+1) = (F(k,i,j) - a*u_z(i,j,k))*adi
        u_z(i,j,k+1) = u_z(i,j,k+1)/dz_deta(k+1)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


       ! CALL COPY(u_z,u_cur)

        END SUBROUTINE
        !______________________________________________________________________

















        !______________________________________________________________________
        SUBROUTINE TURKEL_FORWARD_Z(u_cur,u_z)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NZSB,NZSF
        USE GRID_SIZE_TIME_STEP, ONLY : dz1,azo,dz,dzi,d_eta_z,d_eta_z_inv
        USE HIXON_FORWARD_PARAMETERS
        USE METRICS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),   INTENT(INOUT) :: u_z
        REAL (KIND=8), DIMENSION(NZB:NZF-1,NXL:NXR,NYB:NYT), SAVE :: F

        REAL (KIND=8) :: c
        INTEGER :: i,j,k

!$OMP PARALLEL SHARED(u_cur,F) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        do k = NZB,NZF-1
        F(k,i,j) = (u_cur(i,j,k+1) - u_cur(i,j,k))*d_eta_z_inv
        !dz = dz/azo
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

      
!$OMP PARALLEL SHARED(u_cur,u_z) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        u_z(i,j,NZF) = -4.0d0*u_cur(i,j,NZF-1) + u_cur(i,j,NZF-2) +3.0d0*u_cur(i,j,NZF)
        u_z(i,j,NZF) = 0.5d0 * u_z(i,j,NZF)*d_eta_z_inv
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL



!$OMP PARALLEL SHARED(F,u_z) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        do k = NZF-1,NZB,-1
        u_z(i,j,k) = (F(k,i,j) - a*u_z(i,j,k+1))*adi
        u_z(i,j,k) = u_z(i,j,k)/dz_deta(k)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        !CALL COPY(u_z,u_cur)

        END SUBROUTINE
        !______________________________________________________________________





