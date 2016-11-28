

        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 01-02-2011        Last Modified on 01-02-2011

        !_____________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE TEMPERATURE1

	USE FLOW_PARAMETERS, ONLY : R
        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZF,NZB
	USE continuity,      ONLY : rho_int
        USE energy,          ONLY : T
        USE pressure,        ONLY : p_int
  

        IMPLICIT NONE

        INTEGER :: i,j,k


        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        T(i,j,k) = p_int(i,j,k)/(R*rho_int(i,j,k))

        end do
        end do
        end do

  

        END SUBROUTINE
        !______________________________________________________________________


        !______________________________________________________________________
        SUBROUTINE TEMPERATURE2

	USE FLOW_PARAMETERS, ONLY : R
        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZF,NZB
	USE continuity,      ONLY : rho_new
        USE energy,          ONLY : T
        USE pressure,        ONLY : p_new
  

        IMPLICIT NONE

        INTEGER :: i,j,k
   

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        T(i,j,k) = p_new(i,j,k)/(R*rho_new(i,j,k))

        end do
        end do
        end do

  

        END SUBROUTINE
        !______________________________________________________________________




