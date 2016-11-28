        !______________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 11-07-2009        Last Modified on 10-05-2010

        !______________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE TOTAL_VALUE

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity     , ONLY : rho_cur,rho_tot
        USE x_momentum     , ONLY : u_cur, u_tot,rhou_tot
        USE y_momentum     , ONLY : v_cur, v_tot,rhov_tot
        USE z_momentum     , ONLY : w_cur, w_tot,rhow_tot
        USE pressure       , ONLY : p_cur, p_tot
        USE energy         , ONLY : e_cur, e_tot
        USE cross_momentum , ONLY : uv_tot,uw_tot,vw_tot,turb_ke_tot
        USE lighthill_source


        INTEGER :: i,j,k

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        rho_tot(i,j,k)   = rho_tot(i,j,k) + rho_cur(i,j,k)
        u_tot(i,j,k)     = u_tot(i,j,k) + u_cur(i,j,k)
        v_tot(i,j,k)     = v_tot(i,j,k) + v_cur(i,j,k)
        w_tot(i,j,k)     = w_tot(i,j,k) + w_cur(i,j,k)
        p_tot(i,j,k)     = p_tot(i,j,k) + p_cur(i,j,k)
        e_tot(i,j,k)     = e_tot(i,j,k) + e_cur(i,j,k)
        lh_s_tot(i,j,k)  = lh_s_tot(i,j,k) + lh_s(i,j,k)

        uv_tot(i,j,k)    = uv_tot(i,j,k) + u_cur(i,j,k)*v_cur(i,j,k)
        uw_tot(i,j,k)    = uw_tot(i,j,k) + u_cur(i,j,k)*w_cur(i,j,k)
        vw_tot(i,j,k)    = vw_tot(i,j,k) + v_cur(i,j,k)*w_cur(i,j,k)


        rhou_tot(i,j,k)     = rhou_tot(i,j,k) + rho_cur(i,j,k)*u_cur(i,j,k)
        rhov_tot(i,j,k)     = rhov_tot(i,j,k) + rho_cur(i,j,k)*v_cur(i,j,k)
        rhow_tot(i,j,k)     = rhow_tot(i,j,k) + rho_cur(i,j,k)*w_cur(i,j,k)
        
        turb_ke_tot(i,j,k) = turb_ke_tot(i,j,k) + rho_cur(i,j,k) * (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2)

        end do
        end do
        end do


        END SUBROUTINE
        !______________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE MEAN_VALUE(Nu)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity     , ONLY : rho_tot !,rho_mean 
        USE x_momentum     , ONLY : u_tot,rhou_tot   !, u_mean
        USE y_momentum     , ONLY : v_tot,rhov_tot   !, v_mean
        USE z_momentum     , ONLY : w_tot,rhow_tot   !, w_mean
        USE pressure       , ONLY : p_tot   !, p_mean
        USE energy         , ONLY : e_tot   !, e_mean
        USE cross_momentum , ONLY : uv_tot,uw_tot,vw_tot,turb_ke_tot
        USE lighthill_source
        !USE vorticity   , ONLY : omega_tot, omega_mean

        INTEGER, INTENT(IN) :: Nu
        INTEGER :: i,j,k

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        !rho_mean(i,j,k)   = rho_tot(i,j,k)/Nu
        !u_mean(i,j,k)     = u_tot(i,j,k)/Nu
        !v_mean(i,j,k)     = v_tot(i,j,k)/Nu
        !w_mean(i,j,k)     = w_tot(i,j,k)/Nu
        !p_mean(i,j,k)     = p_tot(i,j,k)/Nu
        !e_mean(i,j,k)     = e_tot(i,j,k)/Nu


        rho_tot(i,j,k)   = rho_tot(i,j,k)/Nu
        u_tot(i,j,k)     = u_tot(i,j,k)/Nu
        v_tot(i,j,k)     = v_tot(i,j,k)/Nu
        w_tot(i,j,k)     = w_tot(i,j,k)/Nu
        p_tot(i,j,k)     = p_tot(i,j,k)/Nu
        e_tot(i,j,k)     = e_tot(i,j,k)/Nu
        lh_s_tot(i,j,k)  = lh_s_tot(i,j,k)/Nu


        uv_tot(i,j,k)    = uv_tot(i,j,k)/Nu
        uw_tot(i,j,k)    = uw_tot(i,j,k)/Nu
        vw_tot(i,j,k)    = vw_tot(i,j,k)/Nu

        rhou_tot(i,j,k)     = rhou_tot(i,j,k)/Nu
        rhov_tot(i,j,k)     = rhov_tot(i,j,k)/Nu
        rhow_tot(i,j,k)     = rhow_tot(i,j,k)/Nu
 
        turb_ke_tot(i,j,k) = turb_ke_tot(i,j,k)/Nu

        end do
        end do
        end do


        END SUBROUTINE
        !______________________________________________________________________


