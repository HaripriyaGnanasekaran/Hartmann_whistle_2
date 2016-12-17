

        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 26-08-2009        Last Modified on 10-05-2010

        !_____________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE VORTICITY1

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE VORTICITY,       ONLY : omegax,omegay,omegaz
        USE x_momentum,      ONLY : u_y,u_z
        USE y_momentum,      ONLY : v_x,v_z
        USE z_momentum,      ONLY : w_x,w_y

        IMPLICIT NONE

        INTEGER :: i,j,k
        !REAL (KIND=8) :: vort

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        omegax(i,j,k) = v_z(i,j,k) - w_y(i,j,k)
        omegay(i,j,k) = w_x(i,j,k) - u_z(i,j,k)
        omegaz(i,j,k) = u_y(i,j,k) - v_x(i,j,k)

        end do
        end do
        end do

        !k = (NZB+NZF)/2
        !do i = NXL,NXR
        !do j = NYB,NYT
        !omega_xyprof(i,j) = v_x(i,j,k) - u_y(i,j,k)
        !end do
        !end do

        END SUBROUTINE
        !______________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE VORTICITY2

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE VORTICITY,       ONLY : omegax,omegay,omegaz
        USE x_momentum,      ONLY : u_y,u_z
        USE y_momentum,      ONLY : v_x,v_z
        USE z_momentum,      ONLY : w_x,w_y

        IMPLICIT NONE

        INTEGER :: i,j,k
        !REAL(KIND = 8) :: vort

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        omegax(i,j,k) = (omegax(i,j,k) + v_z(i,j,k) - w_y(i,j,k)) * 0.5d0
        omegay(i,j,k) = (omegay(i,j,k) + w_x(i,j,k) - u_z(i,j,k)) * 0.5d0
        omegaz(i,j,k) = (omegaz(i,j,k) + u_y(i,j,k) - v_x(i,j,k)) * 0.5d0

        end do
        end do
        end do

        !k = (NZB+NZF)/2
        !do i = NXL,NXR
        !do j = NYB,NYT
        !omega_xyprof(i,j) = (omega_xyprof(i,j) + (v_x(i,j,k) - u_y(i,j,k)))*0.5d0
        !end do
        !end do

        END SUBROUTINE
        !_______________________________________________________________________


