        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 09-11-2009        Last Modified on 10-05-2010

        !_____________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE SPL_VALUE(l)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE pressure    , ONLY : p_tot,p_ts !,p_ms, p_mean
        USE sound_field , ONLY : spl_f

        IMPLICIT NONE
        
        INTEGER, INTENT(IN) :: l
        INTEGER :: i,j,k
        REAL(KIND=8) :: p_val
        REAL(KIND=8) :: ps_val
        REAL(KIND=8),PARAMETER :: i_ref = 10.0d0**(-12)
        REAL(KIND=8),PARAMETER :: rho_ref = 1.25d0
        REAL(KIND=8),PARAMETER :: c_ref = 340.0d0
        REAL(KIND=8),PARAMETER :: p_ref = 1823.85d0/101325.0D0 * 20.0d0*10.0d0**(-6)

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        p_val  = p_tot(i,j,k)/l
        ps_val = p_ts(i,j,k)/l
        ps_val = abs(ps_val - p_val**2)
        !ps_val = ps_val * 0.5d0 /(rho_ref*c_ref)

        !spl_f(i,j,k) = 10.0d0 * log10(ps_val/i_ref)
        ps_val = sqrt(ps_val)
        spl_f(i,j,k) = 20.0d0 * log10(ps_val/p_ref)

        end do
        end do
        end do


        END SUBROUTINE
        !______________________________________________________________________


