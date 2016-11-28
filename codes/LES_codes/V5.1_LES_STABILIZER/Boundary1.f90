        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 11-07-2009        Last Modified on 01-07-2010

        !_____________________________________________________________________




        !_____________________________________________________________________
        SUBROUTINE BOUNDARY1(ts)
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE FLOW_PARAMETERS,     ONLY : rho_atm,p_atm,v_jet        
        USE FLOW_PARAMETERS,     ONLY : c,gama,mu
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dy1,dz1,dt

        USE continuity, ONLY : rho_cur,drho1
        USE continuity, ONLY : rho_x,rho_y,rho_z

        USE x_momentum, ONLY : rhou_cur,drhou1
        USE x_momentum, ONLY : rhou_x!,rhou_y,rhou_z
        USE x_momentum, ONLY : vrhou_xx,vrhou_yy,vrhou_zz
        USE x_momentum, ONLY : rhouu_x,rhouv_y,rhouw_z
        USE x_momentum, ONLY : u_cur,u_x,u_y,u_z

        USE y_momentum, ONLY : rhov_cur,drhov1
        USE y_momentum, ONLY : rhov_y!,rhov_x,rhov_z
        USE y_momentum, ONLY : vrhov_xx,vrhov_yy,vrhov_zz
        USE y_momentum, ONLY : rhovu_x,rhovv_y,rhovw_z
        USE y_momentum, ONLY : v_cur,v_x,v_y,v_z

        USE z_momentum, ONLY : rhow_cur,drhow1
        USE z_momentum, ONLY : rhow_z!,rhow_y,rhow_x
        USE z_momentum, ONLY : vrhow_xx,vrhow_yy,vrhow_zz
        USE z_momentum, ONLY : rhowu_x,rhowv_y,rhoww_z
        USE z_momentum, ONLY : w_cur,w_x,w_y,w_z

        USE pressure,   ONLY : p_cur,p_x,p_y,p_z

        USE energy,     ONLY : e_cur,de1
        USE energy,     ONLY : vke_xx,vke_yy,vke_zz,efl_x,efl_y,efl_z

        IMPLICIT NONE
        
        INTEGER, INTENT(INOUT) :: ts
        INTEGER :: i,j,k

        REAL(KIND=8) :: L1,L2,L3,L4,L5,x,expx
        REAL(KIND=8) :: dp,du,dv,dw,v_sqr,intr

        L1      = 0.0d0
        L2      = 0.0d0
        L3      = 0.0d0
        L4      = 0.0d0
        L5      = 0.0d0
        x       = 0.0d0
        expx    = 0.0d0
    



        !______________________________________________________________________
        !FACE BOUNDARIES

        
        !NYB BOTTOM FACE BOUNDARY

        j = NYB

        do k = NZB+1,NZF-1
        do i = NXL+1,NXR-1

        c = sqrt(gama*p_cur(i,j,k)/rho_cur(i,j,k))

            
        if ((v_cur(i,j,k)-c) >= 0.0d0) then
        L1 = 0.0d0
        else
        L1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_cur(i,j,k)*c*v_y(i,j,k))
        end if

        if(v_cur(i,j,k) >= 0.0d0) then
        L2 = 0.0d0
        L3 = 0.0d0
        L4 = 0.0d0
        else
        L2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        L3 = v_cur(i,j,k)*u_y(i,j,k)
        L4 = v_cur(i,j,k)*w_y(i,j,k)
        end if

        if((v_cur(i,j,k)+c)>= 0.0d0) then
        L5 = 0.0d0
        else
        L5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_cur(i,j,k)*c*v_y(i,j,k))
        end if

        drho1(i,j,k) = (L2+0.5d0*(L5+L1))/(c*c) + u_cur(i,j,k)*rho_x(i,j,k) 
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*u_x(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + w_cur(i,j,k)*rho_z(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*w_z(i,j,k)

        dv = 0.5d0*(L5-L1)/(rho_cur(i,j,k)*c) + u_cur(i,j,k)*v_x(i,j,k)
        dv = dv + w_cur(i,j,k)*v_z(i,j,k)

        du = L3 + p_x(i,j,k)/rho_cur(i,j,k) + u_cur(i,j,k)*u_x(i,j,k)
        du = du + w_cur(i,j,k)*u_z(i,j,k)

        dw = L4 + p_z(i,j,k)/rho_cur(i,j,k) + u_cur(i,j,k)*w_x(i,j,k)
        dw = dw + w_cur(i,j,k)*w_z(i,j,k)

        dp = (L5 + L1)*0.5d0 + u_cur(i,j,k)*p_x(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*u_x(i,j,k)
        dp = dp + w_cur(i,j,k)*p_z(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*w_z(i,j,k)


        drhou1(i,j,k) = u_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*du

        drhov1(i,j,k) = v_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dv

        drhow1(i,j,k) = w_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du +  v_cur(i,j,k)*dv + w_cur(i,j,k)*dw) *rho_cur(i,j,k)
        de1(i,j,k) = v_sqr*drho1(i,j,k) + intr + dp/(gama-1.0d0)

        end do
        end do




        !NYT TOP FACE BOUNDARY

        j = NYT

        do k = NZB+1,NZF-1
        do i = NXL+1,NXR-1

        c = sqrt(gama*p_cur(i,j,k)/rho_cur(i,j,k))

            
        if ((v_cur(i,j,k)-c) <= 0.0d0) then
        L1 = 0.0d0
        else
        L1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_cur(i,j,k)*c*v_y(i,j,k))
        end if

        if(v_cur(i,j,k) <= 0.0d0) then
        L2 = 0.0d0
        L3 = 0.0d0
        L4 = 0.0d0
        else
        L2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        L3 = v_cur(i,j,k)*u_y(i,j,k)
        L4 = v_cur(i,j,k)*w_y(i,j,k)
        end if

        if((v_cur(i,j,k)+c)<= 0.0d0) then
        L5 = 0.0d0
        else
        L5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_cur(i,j,k)*c*v_y(i,j,k))
        end if

        drho1(i,j,k) = (L2+0.5d0*(L5+L1))/(c*c) + u_cur(i,j,k)*rho_x(i,j,k) 
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*u_x(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + w_cur(i,j,k)*rho_z(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*w_z(i,j,k)

        dv = 0.5d0*(L5-L1)/(rho_cur(i,j,k)*c) + u_cur(i,j,k)*v_x(i,j,k)
        dv = dv + w_cur(i,j,k)*v_z(i,j,k)

        du = L3 + p_x(i,j,k)/rho_cur(i,j,k) + u_cur(i,j,k)*u_x(i,j,k)
        du = du + w_cur(i,j,k)*u_z(i,j,k)

        dw = L4 + p_z(i,j,k)/rho_cur(i,j,k) + u_cur(i,j,k)*w_x(i,j,k)
        dw = dw + w_cur(i,j,k)*w_z(i,j,k)

        dp = (L5 + L1)*0.5d0 + u_cur(i,j,k)*p_x(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*u_x(i,j,k)
        dp = dp + w_cur(i,j,k)*p_z(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*w_z(i,j,k)


        drhou1(i,j,k) = u_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*du

        drhov1(i,j,k) = v_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dv

        drhow1(i,j,k) = w_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du +  v_cur(i,j,k)*dv + w_cur(i,j,k)*dw) *rho_cur(i,j,k)
        de1(i,j,k) = v_sqr*drho1(i,j,k) + intr + dp/(gama-1.0d0)

        end do
        end do


        !NXL LEFT FACE BOUNDARY (INFLOW)

        i = NXL

        do k = NZB+1,NZF-1
        do j = NYB+1,NYT-1

        c = sqrt(gama*p_cur(i,j,k)/rho_cur(i,j,k))

 
        if ((u_cur(i,j,k) - c) >= 0.0d0) then
        L1 = 0.0d0
        else
        L1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_cur(i,j,k)*c*u_x(i,j,k))
        end if

        if (u_cur(i,j,k) >= 0.0d0) then
        L2 = 0.0d0
        L3 = 0.0d0
        L4 = 0.0d0
        else
        L2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        L3 = u_cur(i,j,k)*v_x(i,j,k)
        L4 = u_cur(i,j,k)*w_x(i,j,k)
        end if

        if ((u_cur(i,j,k) + c) >= 0.0d0) then
        L5 = 0.0d0
        else
        L5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_cur(i,j,k)*c*u_x(i,j,k))
        end if
        
        drho1(i,j,k) = (L2 + 0.5d0*(L1+L5))/(c*c) + v_cur(i,j,k)*rho_y(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*v_y(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + w_cur(i,j,k)*rho_z(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*w_z(i,j,k)

        du = 0.5d0*(L5-L1)/(rho_cur(i,j,k)*c) + v_cur(i,j,k)*u_y(i,j,k)
        du = du + w_cur(i,j,k)*u_z(i,j,k)

        dv = L3 + p_y(i,j,k)/rho_cur(i,j,k)+ v_cur(i,j,k)*v_y(i,j,k)
        dv = dv + w_cur(i,j,k)*v_z(i,j,k)

        dw = L4 + p_z(i,j,k)/rho_cur(i,j,k)+ v_cur(i,j,k)*w_y(i,j,k)
        dw = dw + w_cur(i,j,k)*w_z(i,j,k)

        dp = 0.5d0*(L5 + L1) + v_cur(i,j,k)*p_y(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*v_y(i,j,k)
        dp = dp + w_cur(i,j,k)*p_z(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*w_z(i,j,k)


        drhou1(i,j,k) = u_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*du

        drhov1(i,j,k) = v_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dv

        drhow1(i,j,k) = w_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_cur(i,j,k)
        de1(i,j,k) = v_sqr*drho1(i,j,k) + intr + dp/(gama-1.0d0)

        end do
        end do





        !NXR RIGHT FACE BOUNDARY (OUTFLOW)

        i = NXR

        do k = NZB+1,NZF-1
        do j = NYB+1,NYT-1

        c = sqrt(gama*p_cur(i,j,k)/rho_cur(i,j,k))

 
        if ((u_cur(i,j,k) - c) <= 0.0d0) then
        L1 = 0.0d0
        else
        L1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_cur(i,j,k)*c*u_x(i,j,k))
        end if

        if (u_cur(i,j,k) <= 0.0d0) then
        L2 = 0.0d0
        L3 = 0.0d0
        L4 = 0.0d0
        else
        L2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        L3 = u_cur(i,j,k)*v_x(i,j,k)
        L4 = u_cur(i,j,k)*w_x(i,j,k)
        end if

        if ((u_cur(i,j,k) + c) <= 0.0d0) then
        L5 = 0.0d0
        else
        L5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_cur(i,j,k)*c*u_x(i,j,k))
        end if
        
        drho1(i,j,k) = (L2 + 0.5d0*(L1+L5))/(c*c) + v_cur(i,j,k)*rho_y(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*v_y(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + w_cur(i,j,k)*rho_z(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*w_z(i,j,k)

        du = 0.5d0*(L5-L1)/(rho_cur(i,j,k)*c) + v_cur(i,j,k)*u_y(i,j,k)
        du = du + w_cur(i,j,k)*u_z(i,j,k)

        dv = L3 + p_y(i,j,k)/rho_cur(i,j,k) + v_cur(i,j,k)*v_y(i,j,k)
        dv = dv + w_cur(i,j,k)*v_z(i,j,k)

        dw = L4 + p_z(i,j,k)/rho_cur(i,j,k) + v_cur(i,j,k)*w_y(i,j,k)
        dw = dw + w_cur(i,j,k)*w_z(i,j,k)

        dp = 0.5d0*(L5 + L1) + v_cur(i,j,k)*p_y(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*v_y(i,j,k)
        dp = dp + w_cur(i,j,k)*p_z(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*w_z(i,j,k)


        drhou1(i,j,k) = u_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*du

        drhov1(i,j,k) = v_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dv

        drhow1(i,j,k) = w_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_cur(i,j,k)
        de1(i,j,k) = v_sqr*drho1(i,j,k) + intr + dp/(gama-1.0d0)

        end do
        end do

       



        !NZB BACK FACE BOUNDARY

        k = NZB

        do j = NYB+1,NYT-1
        do i = NXL+1,NXR-1

        c = sqrt(gama*p_cur(i,j,k)/rho_cur(i,j,k))

            
        if ((w_cur(i,j,k)-c) >= 0.0d0) then
        L1 = 0.0d0
        else
        L1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_cur(i,j,k)*c*w_z(i,j,k))
        end if

        if(w_cur(i,j,k) >= 0.0d0) then
        L2 = 0.0d0
        L3 = 0.0d0
        L4 = 0.0d0
        else
        L2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        L3 = w_cur(i,j,k)*u_z(i,j,k)
        L4 = w_cur(i,j,k)*v_z(i,j,k)
        end if

        if((w_cur(i,j,k)+c)>= 0.0d0) then
        L5 = 0.0d0
        else
        L5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_cur(i,j,k)*c*w_z(i,j,k))
        end if

        drho1(i,j,k) = (L2+0.5d0*(L5+L1))/(c*c) + u_cur(i,j,k)*rho_x(i,j,k) 
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*u_x(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + v_cur(i,j,k)*rho_y(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*v_y(i,j,k)

        dw = 0.5d0*(L5-L1)/(rho_cur(i,j,k)*c) + u_cur(i,j,k)*w_x(i,j,k)
        dw = dw + v_cur(i,j,k)*w_y(i,j,k)

        du = L3 + p_x(i,j,k)/rho_cur(i,j,k)+ u_cur(i,j,k)*u_x(i,j,k)
        du = du + v_cur(i,j,k)*u_y(i,j,k)

        dv = L4 + p_y(i,j,k)/rho_cur(i,j,k)+ u_cur(i,j,k)*v_x(i,j,k)
        dv = dv + v_cur(i,j,k)*v_y(i,j,k)

        dp = (L5 + L1)*0.5d0 + u_cur(i,j,k)*p_x(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*u_x(i,j,k)
        dp = dp + v_cur(i,j,k)*p_y(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*v_y(i,j,k)


        drhou1(i,j,k) = u_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*du

        drhov1(i,j,k) = v_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dv

        drhow1(i,j,k) = w_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du +  v_cur(i,j,k)*dv + w_cur(i,j,k)*dw) *rho_cur(i,j,k)
        de1(i,j,k) = v_sqr*drho1(i,j,k) + intr + dp/(gama-1.0d0)

        end do
        end do



        !NZF FRONT FACE BOUNDARY

        k = NZF

        do j = NYB+1,NYT-1
        do i = NXL+1,NXR-1

        c = sqrt(gama*p_cur(i,j,k)/rho_cur(i,j,k))

            
        if ((w_cur(i,j,k)-c) <= 0.0d0) then
        L1 = 0.0d0
        else
        L1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_cur(i,j,k)*c*w_z(i,j,k))
        end if

        if(w_cur(i,j,k) <= 0.0d0) then
        L2 = 0.0d0
        L3 = 0.0d0
        L4 = 0.0d0
        else
        L2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        L3 = w_cur(i,j,k)*u_z(i,j,k)
        L4 = w_cur(i,j,k)*v_z(i,j,k)
        end if

        if((w_cur(i,j,k)+c)<= 0.0d0) then
        L5 = 0.0d0
        else
        L5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_cur(i,j,k)*c*w_z(i,j,k))
        end if

        drho1(i,j,k) = (L2+0.5d0*(L5+L1))/(c*c) + u_cur(i,j,k)*rho_x(i,j,k) 
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*u_x(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + v_cur(i,j,k)*rho_y(i,j,k)
        drho1(i,j,k) = drho1(i,j,k) + rho_cur(i,j,k)*v_y(i,j,k)

        dw = 0.5d0*(L5-L1)/(rho_cur(i,j,k)*c) + u_cur(i,j,k)*w_x(i,j,k)
        dw = dw + v_cur(i,j,k)*w_y(i,j,k)

        du = L3 + p_x(i,j,k)/rho_cur(i,j,k)+ u_cur(i,j,k)*u_x(i,j,k)
        du = du + v_cur(i,j,k)*u_y(i,j,k)

        dv = L4 + p_y(i,j,k)/rho_cur(i,j,k)+ u_cur(i,j,k)*v_x(i,j,k)
        dv = dv + v_cur(i,j,k)*v_y(i,j,k)

        dp = (L5 + L1)*0.5d0 + u_cur(i,j,k)*p_x(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*u_x(i,j,k)
        dp = dp + v_cur(i,j,k)*p_y(i,j,k)
        dp = dp + gama*p_cur(i,j,k)*v_y(i,j,k)


        drhou1(i,j,k) = u_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*du

        drhov1(i,j,k) = v_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dv

        drhow1(i,j,k) = w_cur(i,j,k)*drho1(i,j,k) + rho_cur(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du +  v_cur(i,j,k)*dv + w_cur(i,j,k)*dw) *rho_cur(i,j,k)
        de1(i,j,k) = v_sqr*drho1(i,j,k) + intr + dp/(gama-1.0d0)

        end do
        end do


        END SUBROUTINE
        !______________________________________________________________________




