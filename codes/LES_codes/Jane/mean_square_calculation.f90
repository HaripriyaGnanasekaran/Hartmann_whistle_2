        !______________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 11-07-2009        Last Modified on 10-05-2010

        !______________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE TOTAL_SQUARE_VALUE

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity  , ONLY : rho_cur,rho_ts
        USE x_momentum  , ONLY : u_cur, u_ts
        USE y_momentum  , ONLY : v_cur, v_ts
        USE z_momentum  , ONLY : w_cur, w_ts
        USE pressure    , ONLY : p_cur, p_ts
        USE lighthill_source

        INTEGER :: i,j,k

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        rho_ts(i,j,k)   = rho_ts(i,j,k) + rho_cur(i,j,k)**2
        u_ts(i,j,k)     = u_ts(i,j,k) + u_cur(i,j,k)**2
        v_ts(i,j,k)     = v_ts(i,j,k) + v_cur(i,j,k)**2
        w_ts(i,j,k)     = w_ts(i,j,k) + w_cur(i,j,k)**2
        p_ts(i,j,k)     = p_ts(i,j,k) + p_cur(i,j,k)**2
        lh_s_ts(i,j,k)  = lh_s_ts(i,j,k) + lh_s(i,j,k)**2
        
        end do
        end do
        end do


        END SUBROUTINE
        !______________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE MEAN_SQUARE_VALUE(Nu)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity  , ONLY : rho_ts !,rho_ms
        USE x_momentum  , ONLY : u_ts   !, u_ms
        USE y_momentum  , ONLY : v_ts   !,v_ms
        USE z_momentum  , ONLY : w_ts   !,w_ms
        USE pressure    , ONLY : p_ts   !,p_ms
        USE lighthill_source

        INTEGER, INTENT(IN) :: Nu
        INTEGER :: i,j,k

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        rho_ts(i,j,k)   = rho_ts(i,j,k)/Nu
        u_ts(i,j,k)     = u_ts(i,j,k)/Nu
        v_ts(i,j,k)     = v_ts(i,j,k)/Nu
        w_ts(i,j,k)     = w_ts(i,j,k)/Nu
        p_ts(i,j,k)     = p_ts(i,j,k)/Nu
        lh_s_ts(i,j,k)  = lh_s_ts(i,j,k)/Nu
        
        end do
        end do
        end do


        END SUBROUTINE
        !______________________________________________________________________


