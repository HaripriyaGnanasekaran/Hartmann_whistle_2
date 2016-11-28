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
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dx,axo,dxi
        USE HIXON_BACKWARD_PARAMETERS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)    :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_x
        REAL (KIND=8), DIMENSION(NXL:NXR-1,NYB:NYT,NZB:NZF), SAVE   :: F
        REAL (KIND=8) :: c
        INTEGER :: i,j,k


!$OMP PARALLEL SHARED(u_cur,u_x,F,dxi) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR-1
        F(i,j,k) = (u_cur(i+1,j,k) - u_cur(i,j,k))*dxi(i)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

        !u_x(NXL,j,k) = b1*u_cur(NXL,j,k)  +  b2*u_cur(NXL+1,j,k)  +  b3*u_cur(NXL+2,j,k)
        !u_x(NXL,j,k) = u_x(NXL,j,k) + b4*u_cur(NXL+3,j,k) + b5*u_cur(NXL+5,j,k)
        !u_x(NXL,j,k) = u_x(NXL,j,k) / dx1

        dx = dx1

        !u_x(NXL,j,k) = 0.0d0 - u_cur(NXL,j,k)*(1.0d0 - c*c) - c*c*u_cur(NXL+2,j,k)
        !u_x(NXL,j,k) = u_x(NXL,j,k) + u_cur(NXL+1,j,k)
        !u_x(NXL,j,k) = u_x(NXL,j,k) / (dx*(1.0d0-c))

!$OMP PARALLEL SHARED(u_cur,u_x) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        u_x(NXL,j,k) = 4.0d0*u_cur(NXL+1,j,k) - u_cur(NXL+2,j,k) - 3.0d0*u_cur(NXL,j,k)
        u_x(NXL,j,k) = 0.5d0 * u_x(NXL,j,k)/dx1
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
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL



        !u_cur = u_x !CALL COPY(u_x,u_cur)

        END SUBROUTINE
        !______________________________________________________________________









        !______________________________________________________________________
        SUBROUTINE TURKEL_FORWARD_X(u_cur,u_x)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NXO,NXB,NYSB,NYST,NZSB,NZSF        
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dx,axo,dxi
        USE HIXON_FORWARD_PARAMETERS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)    :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_x
        REAL (KIND=8), DIMENSION(NXL:NXR-1,NYB:NYT,NZB:NZF), SAVE    :: F

        INTEGER :: i,j,k
        REAL (KIND=8) :: c



!$OMP PARALLEL SHARED(u_cur,F,dxi) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR-1
        F(i,j,k) = (u_cur(i+1,j,k) - u_cur(i,j,k))*dxi(i)
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

        dx = dx1*axo**(NXR-NXO)
        c = axo/(1.0d0 + axo)

!$OMP PARALLEL SHARED(u_cur,u_x) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        u_x(NXR,j,k) = u_cur(NXR,j,k)*(1.0d0 - c*c) + c*c*u_cur(NXR-2,j,k)
        u_x(NXR,j,k) = u_x(NXR,j,k) - u_cur(NXR-1,j,k)
        u_x(NXR,j,k) = u_x(NXR,j,k) / (dx*(1.0d0-c))
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL SHARED(u_x,F) PRIVATE(i)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXR-1,NXL,-1
        u_x(i,j,k) = (F(i,j,k) - a*u_x(i+1,j,k))*adi
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        !u_cur = u_x !CALL COPY(u_x,u_cur)

        END SUBROUTINE
        !______________________________________________________________________









        !______________________________________________________________________
        SUBROUTINE TURKEL_BACKDIF_Y(u_cur,u_y)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NYSB,NYST
        USE GRID_SIZE_TIME_STEP, ONLY : dy1,ayo,dy,dyi
        USE HIXON_BACKWARD_PARAMETERS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)    :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_y
        REAL (KIND=8), DIMENSION(NYB:NYT-1,NXL:NXR,NZB:NZF), SAVE    :: F

        REAL (KIND=8) :: c
        INTEGER :: i,j,k


!$OMP PARALLEL SHARED(u_cur,u_y) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        do j = NYB,NYT-1
        F(j,i,k) = (u_cur(i,j+1,k) - u_cur(i,j,k))*dyi(j)
        !dy = dy/ayo
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        c = ayo/(1.0d0 + ayo)
        dy = dy1 * (ayo**(NYSB-NYB))

!$OMP PARALLEL SHARED(u_cur,u_y) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        u_y(i,NYB,k) = 0.0d0 - u_cur(i,NYB,k)*(1.0d0 - c*c) - c*c*u_cur(i,NYB+2,k)
        u_y(i,NYB,k) = u_y(i,NYB,k) + u_cur(i,NYB+1,k)
        u_y(i,NYB,k) = u_y(i,NYB,k) / (dy*(1.0d0-c))
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
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL



        !u_cur = u_y !CALL COPY(u_y,u_cur)

        END SUBROUTINE
        !______________________________________________________________________

















        !______________________________________________________________________
        SUBROUTINE TURKEL_FORWARD_Y(u_cur,u_y)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NYSB,NYST
        USE GRID_SIZE_TIME_STEP, ONLY : dy1,ayo,dy,dyi
        USE HIXON_FORWARD_PARAMETERS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)    :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_y
        REAL (KIND=8), DIMENSION(NYB:NYT-1,NXL:NXR,NZB:NZF), SAVE    :: F

        REAL (KIND=8) :: c
        INTEGER :: i,j,k


!$OMP PARALLEL SHARED(u_cur,F) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        do j = NYB,NYT-1
        F(j,i,k) = (u_cur(i,j+1,k) - u_cur(i,j,k))*dyi(j)
        !dy = dy/ayo
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

        c = ayo/(1.0d0 + ayo)
        dy = dy1 * (ayo**(NYT-NYST))

!$OMP PARALLEL SHARED(u_cur,u_y) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        u_y(i,NYT,k) = u_cur(i,NYT,k)*(1.0d0 - c*c) + c*c*u_cur(i,NYT-2,k)
        u_y(i,NYT,k) = u_y(i,NYT,k) - u_cur(i,NYT-1,k)
        u_y(i,NYT,k) = u_y(i,NYT,k) / (dy*(1.0d0-c))
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

        j = NYT-1

!$OMP PARALLEL SHARED(F,u_y) PRIVATE(k)
!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do i = NXL,NXR
        do j = NYT-1,NYB,-1
        u_y(i,j,k) = (F(j,i,k) - a*u_y(i,j+1,k))*adi
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        !u_cur = u_y !CALL COPY(u_y,u_cur)

        END SUBROUTINE
        !______________________________________________________________________












        !______________________________________________________________________
        SUBROUTINE TURKEL_BACKDIF_Z(u_cur,u_z)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NZSB,NZSF
        USE GRID_SIZE_TIME_STEP, ONLY : dz1,azo,dz,dzi
        USE HIXON_BACKWARD_PARAMETERS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)    :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_z
        REAL (KIND=8), DIMENSION(NZB:NZF-1,NXL:NXR,NYB:NYT), SAVE    :: F

        REAL (KIND=8) :: c
        INTEGER :: i,j,k


!$OMP PARALLEL SHARED(u_cur,u_z) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        do k = NZB,NZF-1
        F(k,i,j) = (u_cur(i,j,k+1) - u_cur(i,j,k))*dzi(k)
        !dz = dz/azo
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        c = azo/(1.0d0 + azo)
        dz = dz1 * (azo**(NZSB-NZB))

!$OMP PARALLEL SHARED(u_cur,u_z) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        u_z(i,j,NZB) = 0.0d0 - u_cur(i,j,NZB)*(1.0d0 - c*c) - c*c*u_cur(i,j,NZB+2)
        u_z(i,j,NZB) = u_z(i,j,NZB) + u_cur(i,j,NZB+1)
        u_z(i,j,NZB) = u_z(i,j,NZB) / (dz*(1.0d0-c))
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
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL



        !u_cur = u_z !CALL COPY(u_z,u_cur)

        END SUBROUTINE
        !______________________________________________________________________

















        !______________________________________________________________________
        SUBROUTINE TURKEL_FORWARD_Z(u_cur,u_z)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NZSB,NZSF
        USE GRID_SIZE_TIME_STEP, ONLY : dz1,azo,dz,dzi
        USE HIXON_FORWARD_PARAMETERS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(IN)    :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_z
        REAL (KIND=8), DIMENSION(NZB:NZF-1,NXL:NXR,NYB:NYT), SAVE    :: F

        REAL (KIND=8) :: c
        INTEGER :: i,j,k


!$OMP PARALLEL SHARED(u_cur,F) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        do k = NZB,NZF-1
        F(k,i,j) = (u_cur(i,j,k+1) - u_cur(i,j,k))*dzi(k)
        !dz = dz/azo
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL

        c = azo/(1.0d0 + azo)
        dz = dz1 * (azo**(NZF-NZSF))

!$OMP PARALLEL SHARED(u_cur,u_z) PRIVATE(j)
!$OMP DO SCHEDULE(DYNAMIC)
        do j = NYB,NYT
        do i = NXL,NXR
        u_z(i,j,NZF) = u_cur(i,j,NZF)*(1.0d0 - c*c) + c*c*u_cur(i,j,NZF-2)
        u_z(i,j,NZF) = u_z(i,j,NZF) - u_cur(i,j,NZF-1)
        u_z(i,j,NZF) = u_z(i,j,NZF) / (dz*(1.0d0-c))
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
        end do
        end do
        end do
!$OMP END DO
!$OMP END PARALLEL


        !u_cur = u_z !CALL COPY(u_z,u_cur)

        END SUBROUTINE
        !______________________________________________________________________








