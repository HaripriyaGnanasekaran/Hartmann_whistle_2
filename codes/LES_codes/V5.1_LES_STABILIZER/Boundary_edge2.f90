        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 11-07-2009        Last Modified on 27-05-2010

        !_____________________________________________________________________

        !_____________________________________________________________________
        SUBROUTINE BOUNDARY_EDGE2(ch)
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE FLOW_PARAMETERS,     ONLY : rho_atm,p_atm,v_jet
        USE FLOW_PARAMETERS,     ONLY : c,gama
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dy1,dz1,dt

        USE continuity, ONLY : rho_int,drho2
        USE continuity, ONLY : rho_x,rho_y,rho_z

        USE x_momentum, ONLY : rhou_int,drhou2
        USE x_momentum, ONLY : rhou_x!,rhou_y,rhou_z
        USE x_momentum, ONLY : vrhou_xx,vrhou_yy,vrhou_zz
        USE x_momentum, ONLY : rhouu_x,rhouv_y,rhouw_z
        USE x_momentum, ONLY : u_cur,u_x,u_y,u_z

        USE y_momentum, ONLY : rhov_int,drhov2
        USE y_momentum, ONLY : rhov_y!,rhov_x,rhov_z
        USE y_momentum, ONLY : vrhov_xx,vrhov_yy,vrhov_zz
        USE y_momentum, ONLY : rhovu_x,rhovv_y,rhovw_z
        USE y_momentum, ONLY : v_cur,v_x,v_y,v_z

        USE z_momentum, ONLY : rhow_int,drhow2
        USE z_momentum, ONLY : rhow_z!,rhow_y,rhow_x
        USE z_momentum, ONLY : vrhow_xx,vrhow_yy,vrhow_zz
        USE z_momentum, ONLY : rhowu_x,rhowv_y,rhoww_z
        USE z_momentum, ONLY : w_cur,w_x,w_y,w_z

        USE pressure,   ONLY : p_int,p_x,p_y,p_z

        USE energy,     ONLY : e_int,de2
        USE energy,     ONLY : vke_xx,vke_yy,vke_zz,efl_x,efl_y,efl_z

        IMPLICIT NONE
        
        INTEGER, DIMENSION(2), INTENT(IN) :: ch
        INTEGER :: i,j,k

        REAL(KIND=8) :: LX1,LX2,LX3,LX4,LX5,x,expx
        REAL(KIND=8) :: LY1,LY2,LY3,LY4,LY5
        REAL(KIND=8) :: LZ1,LZ2,LZ3,LZ4,LZ5
        REAL(KIND=8) :: dp,du,dv,dw,v_sqr,intr

        LX1      = 0.0d0
        LX2      = 0.0d0
        LX3      = 0.0d0
        LX4      = 0.0d0
        LX5      = 0.0d0
        LY1      = 0.0d0
        LY2      = 0.0d0
        LY3      = 0.0d0
        LY4      = 0.0d0
        LY5      = 0.0d0
        LZ1      = 0.0d0
        LZ2      = 0.0d0
        LZ3      = 0.0d0
        LZ4      = 0.0d0
        LZ5      = 0.0d0
        x       = 0.0d0
        expx    = 0.0d0

        !______________________________________________________________________
        !EDGE TREATMENT

        !NOTE: There are totally 12 edges
        !Edges of faces normal to NXL and NXR make eight
        !Edges along the x direction are four
        !Bottom BACK, Bottom FRONT, Top BACK, Top FRONT

        !LEFT BOTTOM and LEFT TOP EDGES (Axis alond z)
        !Resolved along y direction, x direction derivatives not
        !included. Z direction derivatives included.

        !LEFT BACK and LEFT FRONT (periodic) (Axis along y)
        !Resolved along x direction. y and z direction dervatives
        !included.

        !RIGHT BOTTOM and RIGHT TOP EDGES (Axis alond z)
        !Resolved along y direction, x direction derivatives not
        !included. Z direction derivatives included.

        !RIGHT BACK and RIGHT FRONT (periodic) (Axis along y)
        !Resolved along x direction. y and z direction dervatives
        !included.

        !BOTTOM BACK and BOTTOM FRONT (Axis along x axis)
        !Resolved along y direction, x and z direction derivatives
        !included.

        !TOP BACK and TOP FRONT (Axis along x axis)
        !Resolved along x direction, x and z direction derivatives
        !included.

        !NOTE: Resolved along y,x direction dervatives not included.
        !z direction derivatives included


        !LEFT BOTTOM EDGE
        i = NXL
        j = NYB
        do k = NZB+1,NZF-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((u_cur(i,j,k) - c) >= 0.0d0) then
        LX1 = 0.0d0
        else
        LX1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_int(i,j,k)*c*u_x(i,j,k))
        end if
        if (u_cur(i,j,k) >= 0.0d0) then
        LX2 = 0.0d0
        LX3 = 0.0d0
        LX4 = 0.0d0
        else
        LX2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        LX3 = u_cur(i,j,k)*v_x(i,j,k)
        LX4 = u_cur(i,j,k)*w_x(i,j,k)
        end if
        if ((u_cur(i,j,k) + c) >= 0.0d0) then
        LX5 = 0.0d0
        else
        LX5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_int(i,j,k)*c*u_x(i,j,k))
        end if


        if ((v_cur(i,j,k)-c) >= 0.0d0) then
        LY1 = 0.0d0
        else
        LY1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_int(i,j,k)*c*v_y(i,j,k))
        end if
        if(v_cur(i,j,k) >= 0.0d0) then
        LY2 = 0.0d0
        LY3 = 0.0d0
        LY4 = 0.0d0
        else
        LY2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        LY3 = v_cur(i,j,k)*u_y(i,j,k)
        LY4 = v_cur(i,j,k)*w_y(i,j,k)
        end if
        if((v_cur(i,j,k)+c)>= 0.0d0) then
        LY5 = 0.0d0
        else
        LY5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_int(i,j,k)*c*v_y(i,j,k))
        end if


        drho2(i,j,k) = (LX2+0.5d0*(LX5+LX1))/(c*c) + (LY2+0.5d0*(LY5+LY1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + w_cur(i,j,k)*rho_z(i,j,k)
        drho2(i,j,k) = drho2(i,j,k) + rho_int(i,j,k)*w_z(i,j,k)

        du = 0.5d0*(LX5-LX1)/(rho_int(i,j,k)*c) + LY3 + w_cur(i,j,k)*u_z(i,j,k)

        dv = LX3 + 0.5d0*(LY5-LY1)/(rho_int(i,j,k)*c) + w_cur(i,j,k)*v_z(i,j,k)

        dw = LX4 + LY4 + p_z(i,j,k)/rho_int(i,j,k) + w_cur(i,j,k)*w_z(i,j,k)

        dp = (LX5 + LX1)*0.5d0 + (LY5 + LY1)*0.5d0 
        dp = dp + w_cur(i,j,k)*p_z(i,j,k) + gama*p_int(i,j,k)*w_z(i,j,k)


        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)

        end do




        !LEFT TOP EDGE
        i = NXL
        j = NYT
        do k = NZB+1,NZF-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((u_cur(i,j,k) - c) >= 0.0d0) then
        LX1 = 0.0d0
        else
        LX1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_int(i,j,k)*c*u_x(i,j,k))
        end if
        if (u_cur(i,j,k) >= 0.0d0) then
        LX2 = 0.0d0
        LX3 = 0.0d0
        LX4 = 0.0d0
        else
        LX2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        LX3 = u_cur(i,j,k)*v_x(i,j,k)
        LX4 = u_cur(i,j,k)*w_x(i,j,k)
        end if
        if ((u_cur(i,j,k) + c) >= 0.0d0) then
        LX5 = 0.0d0
        else
        LX5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_int(i,j,k)*c*u_x(i,j,k))
        end if


        if ((v_cur(i,j,k)-c) <= 0.0d0) then
        LY1 = 0.0d0
        else
        LY1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_int(i,j,k)*c*v_y(i,j,k))
        end if
        if(v_cur(i,j,k) <= 0.0d0) then
        LY2 = 0.0d0
        LY3 = 0.0d0
        LY4 = 0.0d0
        else
        LY2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        LY3 = v_cur(i,j,k)*u_y(i,j,k)
        LY4 = v_cur(i,j,k)*w_y(i,j,k)
        end if
        if((v_cur(i,j,k)+c)<= 0.0d0) then
        LY5 = 0.0d0
        else
        LY5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_int(i,j,k)*c*v_y(i,j,k))
        end if


        drho2(i,j,k) = (LX2+0.5d0*(LX5+LX1))/(c*c) + (LY2+0.5d0*(LY5+LY1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + w_cur(i,j,k)*rho_z(i,j,k)
        drho2(i,j,k) = drho2(i,j,k) + rho_int(i,j,k)*w_z(i,j,k)

        du = 0.5d0*(LX5-LX1)/(rho_int(i,j,k)*c) + LY3 + w_cur(i,j,k)*u_z(i,j,k)

        dv = 0.5d0*(LY5-LY1)/(rho_int(i,j,k)*c) + LX3 + w_cur(i,j,k)*v_z(i,j,k)

        dw = LX4 + LY4 + p_z(i,j,k)/rho_int(i,j,k) + w_cur(i,j,k)*w_z(i,j,k)

        dp = (LX5 + LX1)*0.5d0 + (LY5 + LY1)*0.5d0 
        dp = dp + w_cur(i,j,k)*p_z(i,j,k)
        dp = dp + gama*p_int(i,j,k)*w_z(i,j,k)


        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)

        end do




        !LEFT BACK EDGE
        i = NXL
        k = NZB

        do j = NYB+1,NYT-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((u_cur(i,j,k) - c) >= 0.0d0) then
        LX1 = 0.0d0
        else
        LX1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_int(i,j,k)*c*u_x(i,j,k))
        end if
        if (u_cur(i,j,k) >= 0.0d0) then
        LX2 = 0.0d0
        LX3 = 0.0d0
        LX4 = 0.0d0
        else
        LX2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        LX3 = u_cur(i,j,k)*v_x(i,j,k)
        LX4 = u_cur(i,j,k)*w_x(i,j,k)
        end if
        if ((u_cur(i,j,k) + c) >= 0.0d0) then
        LX5 = 0.0d0
        else
        LX5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_int(i,j,k)*c*u_x(i,j,k))
        end if


        if ((w_cur(i,j,k)-c) >= 0.0d0) then
        LZ1 = 0.0d0
        else
        LZ1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_int(i,j,k)*c*w_z(i,j,k))
        end if
        if(w_cur(i,j,k) >= 0.0d0) then
        LZ2 = 0.0d0
        LZ3 = 0.0d0
        LZ4 = 0.0d0
        else
        LZ2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        LZ3 = w_cur(i,j,k)*u_z(i,j,k)
        LZ4 = w_cur(i,j,k)*v_z(i,j,k)
        end if
        if((w_cur(i,j,k)+c)>= 0.0d0) then
        LZ5 = 0.0d0
        else
        LZ5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_int(i,j,k)*c*w_z(i,j,k))
        end if

        drho2(i,j,k) = (LX2+0.5d0*(LX5+LX1))/(c*c) + (LZ2+0.5d0*(LZ5+LZ1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + v_cur(i,j,k)*rho_y(i,j,k)
        drho2(i,j,k) = drho2(i,j,k) + rho_int(i,j,k)*v_y(i,j,k)

        du = 0.5d0*(LX5-LX1)/(rho_int(i,j,k)*c) + LZ3 + v_cur(i,j,k)*u_y(i,j,k)

        dw = 0.5d0*(LZ5-LZ1)/(rho_int(i,j,k)*c) + LX4 + v_cur(i,j,k)*w_y(i,j,k)


        dv = LX3 + LZ4 + p_y(i,j,k)/rho_int(i,j,k) + v_cur(i,j,k)*v_y(i,j,k)

        dp = (LX5 + LX1)*0.5d0 + (LZ5 + LZ1)*0.5d0 
        dp = dp + v_cur(i,j,k)*p_y(i,j,k) + gama*p_int(i,j,k)*v_y(i,j,k)


        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)

        end do




        !LEFT FRONT EDGE
        i = NXL
        k = NZF

        do j = NYB+1,NYT-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((u_cur(i,j,k) - c) >= 0.0d0) then
        LX1 = 0.0d0
        else
        LX1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_int(i,j,k)*c*u_x(i,j,k))
        end if
        if (u_cur(i,j,k) >= 0.0d0) then
        LX2 = 0.0d0
        LX3 = 0.0d0
        LX4 = 0.0d0
        else
        LX2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        LX3 = u_cur(i,j,k)*v_x(i,j,k)
        LX4 = u_cur(i,j,k)*w_x(i,j,k)
        end if
        if ((u_cur(i,j,k) + c) >= 0.0d0) then
        LX5 = 0.0d0
        else
        LX5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_int(i,j,k)*c*u_x(i,j,k))
        end if


        if ((w_cur(i,j,k)-c) <= 0.0d0) then
        LZ1 = 0.0d0
        else
        LZ1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_int(i,j,k)*c*w_z(i,j,k))
        end if
        if(w_cur(i,j,k) <= 0.0d0) then
        LZ2 = 0.0d0
        LZ3 = 0.0d0
        LZ4 = 0.0d0
        else
        LZ2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        LZ3 = w_cur(i,j,k)*u_z(i,j,k)
        LZ4 = w_cur(i,j,k)*v_z(i,j,k)
        end if
        if((w_cur(i,j,k)+c)<= 0.0d0) then
        LZ5 = 0.0d0
        else
        LZ5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_int(i,j,k)*c*w_z(i,j,k))
        end if


        drho2(i,j,k) = (LX2+0.5d0*(LX5+LX1))/(c*c) + (LZ2+0.5d0*(LZ5+LZ1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + v_cur(i,j,k)*rho_y(i,j,k)
        drho2(i,j,k) = drho2(i,j,k) + rho_int(i,j,k)*v_y(i,j,k)

        du = 0.5d0*(LX5-LX1)/(rho_int(i,j,k)*c) + LZ3 + v_cur(i,j,k)*u_y(i,j,k)

        dw = 0.5d0*(LZ5-LZ1)/(rho_int(i,j,k)*c) + LX4 + v_cur(i,j,k)*w_y(i,j,k)

        dv = LX3 + LZ4 + p_y(i,j,k)/rho_int(i,j,k) + v_cur(i,j,k)*v_y(i,j,k)

        dp = (LX5 + LX1)*0.5d0 + (LZ5 + LZ1)*0.5d0 
        dp = dp + v_cur(i,j,k)*p_y(i,j,k) + gama*p_int(i,j,k)*v_y(i,j,k)


        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)


        end do




        !RIGHT BOTTOM EDGE

        i = NXR
        j = NYB

        do k = NZB+1,NZF-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((u_cur(i,j,k) - c) <= 0.0d0) then
        LX1 = 0.0d0
        else
        LX1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_int(i,j,k)*c*u_x(i,j,k))
        end if
        if (u_cur(i,j,k) <= 0.0d0) then
        LX2 = 0.0d0
        LX3 = 0.0d0
        LX4 = 0.0d0
        else
        LX2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        LX3 = u_cur(i,j,k)*v_x(i,j,k)
        LX4 = u_cur(i,j,k)*w_x(i,j,k)
        end if
        if ((u_cur(i,j,k) + c) <= 0.0d0) then
        LX5 = 0.0d0
        else
        LX5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_int(i,j,k)*c*u_x(i,j,k))
        end if


        if ((v_cur(i,j,k)-c) >= 0.0d0) then
        LY1 = 0.0d0
        else
        LY1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_int(i,j,k)*c*v_y(i,j,k))
        end if
        if(v_cur(i,j,k) >= 0.0d0) then
        LY2 = 0.0d0
        LY3 = 0.0d0
        LY4 = 0.0d0
        else
        LY2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        LY3 = v_cur(i,j,k)*u_y(i,j,k)
        LY4 = v_cur(i,j,k)*w_y(i,j,k)
        end if
        if((v_cur(i,j,k)+c)>= 0.0d0) then
        LY5 = 0.0d0
        else
        LY5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_int(i,j,k)*c*v_y(i,j,k))
        end if


        drho2(i,j,k) = (LX2+0.5d0*(LX5+LX1))/(c*c) + (LY2+0.5d0*(LY5+LY1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + w_cur(i,j,k)*rho_z(i,j,k) + rho_int(i,j,k)*w_z(i,j,k)

        dv = 0.5d0*(LY5-LY1)/(rho_int(i,j,k)*c) + LX3 + w_cur(i,j,k)*v_z(i,j,k)

        du = 0.5d0*(LX5-LX1)/(rho_int(i,j,k)*c) + LY3 + w_cur(i,j,k)*u_z(i,j,k)

        dw = LX4 + LY4 + p_z(i,j,k)/rho_int(i,j,k) + w_cur(i,j,k)*w_z(i,j,k)

        dp = (LX5 + LX1)*0.5d0 + (LY5 + LY1)*0.5d0 
        dp = dp + w_cur(i,j,k)*p_z(i,j,k) + gama*p_int(i,j,k)*w_z(i,j,k)


        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)

        end do




        !RIGHT TOP EDGE

        i = NXR
        j = NYT

        do k = NZB+1,NZF-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((u_cur(i,j,k) - c) <= 0.0d0) then
        LX1 = 0.0d0
        else
        LX1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_int(i,j,k)*c*u_x(i,j,k))
        end if
        if (u_cur(i,j,k) <= 0.0d0) then
        LX2 = 0.0d0
        LX3 = 0.0d0
        LX4 = 0.0d0
        else
        LX2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        LX3 = u_cur(i,j,k)*v_x(i,j,k)
        LX4 = u_cur(i,j,k)*w_x(i,j,k)
        end if
        if ((u_cur(i,j,k) + c) <= 0.0d0) then
        LX5 = 0.0d0
        else
        LX5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_int(i,j,k)*c*u_x(i,j,k))
        end if


        if ((v_cur(i,j,k)-c) <= 0.0d0) then
        LY1 = 0.0d0
        else
        LY1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_int(i,j,k)*c*v_y(i,j,k))
        end if
        if(v_cur(i,j,k) <= 0.0d0) then
        LY2 = 0.0d0
        LY3 = 0.0d0
        LY4 = 0.0d0
        else
        LY2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        LY3 = v_cur(i,j,k)*u_y(i,j,k)
        LY4 = v_cur(i,j,k)*w_y(i,j,k)
        end if
        if((v_cur(i,j,k)+c)<= 0.0d0) then
        LY5 = 0.0d0
        else
        LY5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_int(i,j,k)*c*v_y(i,j,k))
        end if


        drho2(i,j,k) = (LX2+0.5d0*(LX5+LX1))/(c*c) + (LY2+0.5d0*(LY5+LY1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + w_cur(i,j,k)*rho_z(i,j,k) + rho_int(i,j,k)*w_z(i,j,k)

        dv = 0.5d0*(LY5-LY1)/(rho_int(i,j,k)*c) + LX3 + w_cur(i,j,k)*v_z(i,j,k)

        du = 0.5d0*(LX5-LX1)/(rho_int(i,j,k)*c) + LY3 + w_cur(i,j,k)*u_z(i,j,k)

        dw = LX4 + LY4 + p_z(i,j,k)/rho_int(i,j,k) + w_cur(i,j,k)*w_z(i,j,k)

        dp = (LX5 + LX1)*0.5d0 + (LY5 + LY1)*0.5d0 
        dp = dp + w_cur(i,j,k)*p_z(i,j,k) + gama*p_int(i,j,k)*w_z(i,j,k)


        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)

        end do


        

        !RIGHT BACK EDGE

        i = NXR
        k = NZB

        do j = NYB+1,NYT-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((u_cur(i,j,k) - c) <= 0.0d0) then
        LX1 = 0.0d0
        else
        LX1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_int(i,j,k)*c*u_x(i,j,k))
        end if
        if (u_cur(i,j,k) <= 0.0d0) then
        LX2 = 0.0d0
        LX3 = 0.0d0
        LX4 = 0.0d0
        else
        LX2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        LX3 = u_cur(i,j,k)*v_x(i,j,k)
        LX4 = u_cur(i,j,k)*w_x(i,j,k)
        end if
        if ((u_cur(i,j,k) + c) <= 0.0d0) then
        LX5 = 0.0d0
        else
        LX5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_int(i,j,k)*c*u_x(i,j,k))
        end if


        if ((w_cur(i,j,k)-c) >= 0.0d0) then
        LZ1 = 0.0d0
        else
        LZ1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_int(i,j,k)*c*w_z(i,j,k))
        end if
        if(w_cur(i,j,k) >= 0.0d0) then
        LZ2 = 0.0d0
        LZ3 = 0.0d0
        LZ4 = 0.0d0
        else
        LZ2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        LZ3 = w_cur(i,j,k)*u_z(i,j,k)
        LZ4 = w_cur(i,j,k)*v_z(i,j,k)
        end if
        if((w_cur(i,j,k)+c)>= 0.0d0) then
        LZ5 = 0.0d0
        else
        LZ5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_int(i,j,k)*c*w_z(i,j,k))
        end if


        drho2(i,j,k) = (LX2+0.5d0*(LX5+LX1))/(c*c) + (LZ2+0.5d0*(LZ5+LZ1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + v_cur(i,j,k)*rho_y(i,j,k) + rho_int(i,j,k)*v_y(i,j,k)

        dw = 0.5d0*(LZ5-LZ1)/(rho_int(i,j,k)*c) + LX4 + v_cur(i,j,k)*w_y(i,j,k)

        du = 0.5d0*(LX5-LX1)/(rho_int(i,j,k)*c) + LZ3 + v_cur(i,j,k)*u_y(i,j,k)

        dv = LX3 + LZ4 + p_y(i,j,k)/rho_int(i,j,k) + v_cur(i,j,k)*v_y(i,j,k)

        dp = (LZ5 + LZ1)*0.5d0 + (LX5 + LX1)*0.5d0 
        dp = dp + v_cur(i,j,k)*p_y(i,j,k) + gama*p_int(i,j,k)*v_y(i,j,k)


        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)

        end do




        !RIGHT FRONT EDGE

        i = NXR
        k = NZF

        do j = NYB+1,NYT-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((u_cur(i,j,k) - c) <= 0.0d0) then
        LX1 = 0.0d0
        else
        LX1 = (u_cur(i,j,k) - c) *(p_x(i,j,k) -rho_int(i,j,k)*c*u_x(i,j,k))
        end if
        if (u_cur(i,j,k) <= 0.0d0) then
        LX2 = 0.0d0
        LX3 = 0.0d0
        LX4 = 0.0d0
        else
        LX2 = u_cur(i,j,k)*(c*c*rho_x(i,j,k) - p_x(i,j,k))
        LX3 = u_cur(i,j,k)*v_x(i,j,k)
        LX4 = u_cur(i,j,k)*w_x(i,j,k)
        end if
        if ((u_cur(i,j,k) + c) <= 0.0d0) then
        LX5 = 0.0d0
        else
        LX5 = (u_cur(i,j,k)+c)*(p_x(i,j,k)+rho_int(i,j,k)*c*u_x(i,j,k))
        end if


        if ((w_cur(i,j,k)-c) <= 0.0d0) then
        LZ1 = 0.0d0
        else
        LZ1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_int(i,j,k)*c*w_z(i,j,k))
        end if
        if(w_cur(i,j,k) <= 0.0d0) then
        LZ2 = 0.0d0
        LZ3 = 0.0d0
        LZ4 = 0.0d0
        else
        LZ2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        LZ3 = w_cur(i,j,k)*u_z(i,j,k)
        LZ4 = w_cur(i,j,k)*v_z(i,j,k)
        end if
        if((w_cur(i,j,k)+c)<= 0.0d0) then
        LZ5 = 0.0d0
        else
        LZ5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_int(i,j,k)*c*w_z(i,j,k))
        end if


        drho2(i,j,k) = (LX2+0.5d0*(LX5+LX1))/(c*c) + (LZ2+0.5d0*(LZ5+LZ1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + v_cur(i,j,k)*rho_y(i,j,k) + rho_int(i,j,k)*v_y(i,j,k)

        dw = 0.5d0*(LZ5-LZ1)/(rho_int(i,j,k)*c) + LX4 + v_cur(i,j,k)*w_y(i,j,k)

        du = 0.5d0*(LX5-LX1)/(rho_int(i,j,k)*c) + LZ3 + v_cur(i,j,k)*u_y(i,j,k)

        dv = LX3 + LZ4 + p_y(i,j,k)/rho_int(i,j,k) + v_cur(i,j,k)*v_y(i,j,k)

        dp = (LX5 + LX1)*0.5d0 + (LZ5 + LZ1)*0.5d0 
        dp = dp + v_cur(i,j,k)*p_y(i,j,k) + gama*p_int(i,j,k)*v_y(i,j,k)


        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)


        end do






        !BOTTOM BACK EDGE

        j = NYB
        k = NZB

        do i = NXL+1,NXR-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))
        
        if ((v_cur(i,j,k)-c) >= 0.0d0) then
        LY1 = 0.0d0
        else
        LY1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_int(i,j,k)*c*v_y(i,j,k))
        end if
        if(v_cur(i,j,k) >= 0.0d0) then
        LY2 = 0.0d0
        LY3 = 0.0d0
        LY4 = 0.0d0
        else
        LY2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        LY3 = v_cur(i,j,k)*u_y(i,j,k)
        LY4 = v_cur(i,j,k)*w_y(i,j,k)
        end if
        if((v_cur(i,j,k)+c)>= 0.0d0) then
        LY5 = 0.0d0
        else
        LY5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_int(i,j,k)*c*v_y(i,j,k))
        end if


        if ((w_cur(i,j,k)-c) >= 0.0d0) then
        LZ1 = 0.0d0
        else
        LZ1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_int(i,j,k)*c*w_z(i,j,k))
        end if
        if(w_cur(i,j,k) >= 0.0d0) then
        LZ2 = 0.0d0
        LZ3 = 0.0d0
        LZ4 = 0.0d0
        else
        LZ2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        LZ3 = w_cur(i,j,k)*u_z(i,j,k)
        LZ4 = w_cur(i,j,k)*v_z(i,j,k)
        end if
        if((w_cur(i,j,k)+c)>= 0.0d0) then
        LZ5 = 0.0d0
        else
        LZ5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_int(i,j,k)*c*w_z(i,j,k))
        end if


        drho2(i,j,k) = (LY2+0.5d0*(LY5+LY1))/(c*c) + (LZ2+0.5d0*(LZ5+LZ1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + rho_int(i,j,k)*u_x(i,j,k) + u_cur(i,j,k)*rho_x(i,j,k)

        dv = 0.5d0*(LY5-LY1)/(rho_int(i,j,k)*c) + LZ4 + u_cur(i,j,k)*v_x(i,j,k)

        du = LY3 + LZ3 + p_x(i,j,k)/rho_int(i,j,k) + u_cur(i,j,k)*u_x(i,j,k)

        dw = 0.5d0*(LZ5-LZ1)/(rho_int(i,j,k)*c) + LY4 + u_cur(i,j,k)*w_x(i,j,k)

        dp = (LY5 + LY1)*0.5d0 + (LZ5 + LZ1)*0.5d0 
        dp = dp + gama*p_int(i,j,k)*u_x(i,j,k) + u_cur(i,j,k)*p_x(i,j,k)


        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)

 
        end do











        !BOTTOM FRONT EDGE
        j = NYB
        k = NZF

        do i = NXL+1,NXR-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))
        
        if ((v_cur(i,j,k)-c) >= 0.0d0) then
        LY1 = 0.0d0
        else
        LY1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_int(i,j,k)*c*v_y(i,j,k))
        end if
        if(v_cur(i,j,k) >= 0.0d0) then
        LY2 = 0.0d0
        LY3 = 0.0d0
        LY4 = 0.0d0
        else
        LY2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        LY3 = v_cur(i,j,k)*u_y(i,j,k)
        LY4 = v_cur(i,j,k)*w_y(i,j,k)
        end if
        if((v_cur(i,j,k)+c)>= 0.0d0) then
        LY5 = 0.0d0
        else
        LY5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_int(i,j,k)*c*v_y(i,j,k))
        end if


        if ((w_cur(i,j,k)-c) <= 0.0d0) then
        LZ1 = 0.0d0
        else
        LZ1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_int(i,j,k)*c*w_z(i,j,k))
        end if
        if(w_cur(i,j,k) <= 0.0d0) then
        LZ2 = 0.0d0
        LZ3 = 0.0d0
        LZ4 = 0.0d0
        else
        LZ2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        LZ3 = w_cur(i,j,k)*u_z(i,j,k)
        LZ4 = w_cur(i,j,k)*v_z(i,j,k)
        end if
        if((w_cur(i,j,k)+c)<= 0.0d0) then
        LZ5 = 0.0d0
        else
        LZ5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_int(i,j,k)*c*w_z(i,j,k))
        end if


        drho2(i,j,k) = (LY2+0.5d0*(LY5+LY1))/(c*c) + (LZ2+0.5d0*(LZ5+LZ1))/(c*c)
        drho2(i,j,k) = drho2(i,j,k) + rho_int(i,j,k)*u_x(i,j,k) + u_cur(i,j,k)*rho_x(i,j,k)

        dv = 0.5d0*(LY5-LY1)/(rho_int(i,j,k)*c) + LZ4 + u_cur(i,j,k)*v_x(i,j,k)

        du = LY3 + LZ3 + p_x(i,j,k)/rho_int(i,j,k) + u_cur(i,j,k)*u_x(i,j,k)

        dw = 0.5d0*(LZ5-LZ1)/(rho_int(i,j,k)*c) + LY4 + u_cur(i,j,k)*w_x(i,j,k)

        dp = (LY5 + LY1)*0.5d0 + (LZ5 + LZ1)*0.5d0 
        dp = dp + gama*p_int(i,j,k)*u_x(i,j,k) + u_cur(i,j,k)*p_x(i,j,k) 

        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)


        end do



        !TOP BACK EDGE
        j = NYT
        k = NZB

        do  i = NXL+1,NXR-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((v_cur(i,j,k)-c) <= 0.0d0) then
        LY1 = 0.0d0
        else
        LY1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_int(i,j,k)*c*v_y(i,j,k))
        end if
        if(v_cur(i,j,k) <= 0.0d0) then
        LY2 = 0.0d0
        LY3 = 0.0d0
        LY4 = 0.0d0
        else
        LY2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        LY3 = v_cur(i,j,k)*u_y(i,j,k)
        LY4 = v_cur(i,j,k)*w_y(i,j,k)
        end if
        if((v_cur(i,j,k)+c)<= 0.0d0) then
        LY5 = 0.0d0
        else
        LY5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_int(i,j,k)*c*v_y(i,j,k))
        end if


        if ((w_cur(i,j,k)-c) >= 0.0d0) then
        LZ1 = 0.0d0
        else
        LZ1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_int(i,j,k)*c*w_z(i,j,k))
        end if
        if(w_cur(i,j,k) >= 0.0d0) then
        LZ2 = 0.0d0
        LZ3 = 0.0d0
        LZ4 = 0.0d0
        else
        LZ2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        LZ3 = w_cur(i,j,k)*u_z(i,j,k)
        LZ4 = w_cur(i,j,k)*v_z(i,j,k)
        end if
        if((w_cur(i,j,k)+c)>= 0.0d0) then
        LZ5 = 0.0d0
        else
        LZ5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_int(i,j,k)*c*w_z(i,j,k))
        end if


        drho2(i,j,k) = (LY2+0.5d0*(LY5+LY1))/(c*c) + (LZ2+0.5d0*(LZ5+LZ1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + rho_int(i,j,k)*u_x(i,j,k) + u_cur(i,j,k)*rho_x(i,j,k)

        dv = 0.5d0*(LY5-LY1)/(rho_int(i,j,k)*c) + LZ4 + u_cur(i,j,k)*v_x(i,j,k)

        du = LY3 + LZ3 + p_x(i,j,k)/rho_int(i,j,k) + u_cur(i,j,k)*u_x(i,j,k)

        dw = 0.5d0*(LZ5-LZ1)/(rho_int(i,j,k)*c) + LY4 + u_cur(i,j,k)*w_x(i,j,k)

        dp = (LY5 + LY1)*0.5d0 + (LZ5 + LZ1)*0.5d0  
        dp = dp + gama*p_int(i,j,k)*u_x(i,j,k) + u_cur(i,j,k)*p_x(i,j,k)

        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)


        end do



        !TOP FRONT EDGE

        j = NYT
        k = NZF

        do i = NXL+1,NXR-1

        c = sqrt(gama*p_int(i,j,k)/rho_int(i,j,k))


        if ((v_cur(i,j,k)-c) <= 0.0d0) then
        LY1 = 0.0d0
        else
        LY1 = (v_cur(i,j,k) - c)*(p_y(i,j,k) - rho_int(i,j,k)*c*v_y(i,j,k))
        end if
        if(v_cur(i,j,k) <= 0.0d0) then
        LY2 = 0.0d0
        LY3 = 0.0d0
        LY4 = 0.0d0
        else
        LY2 = v_cur(i,j,k)*(c*c*rho_y(i,j,k) - p_y(i,j,k))
        LY3 = v_cur(i,j,k)*u_y(i,j,k)
        LY4 = v_cur(i,j,k)*w_y(i,j,k)
        end if
        if((v_cur(i,j,k)+c)<= 0.0d0) then
        LY5 = 0.0d0
        else
        LY5 = (v_cur(i,j,k)+c)*(p_y(i,j,k)+rho_int(i,j,k)*c*v_y(i,j,k))
        end if


        if ((w_cur(i,j,k)-c) <= 0.0d0) then
        LZ1 = 0.0d0
        else
        LZ1 = (w_cur(i,j,k) - c)*(p_z(i,j,k) - rho_int(i,j,k)*c*w_z(i,j,k))
        end if
        if(w_cur(i,j,k) <= 0.0d0) then
        LZ2 = 0.0d0
        LZ3 = 0.0d0
        LZ4 = 0.0d0
        else
        LZ2 = w_cur(i,j,k)*(c*c*rho_z(i,j,k) - p_z(i,j,k))
        LZ3 = w_cur(i,j,k)*u_z(i,j,k)
        LZ4 = w_cur(i,j,k)*v_z(i,j,k)
        end if
        if((w_cur(i,j,k)+c)<= 0.0d0) then
        LZ5 = 0.0d0
        else
        LZ5 = (w_cur(i,j,k)+c)*(p_z(i,j,k) + rho_int(i,j,k)*c*w_z(i,j,k))
        end if


        drho2(i,j,k) = (LY2+0.5d0*(LY5+LY1))/(c*c) + (LZ2+0.5d0*(LZ5+LZ1))/(c*c) 
        drho2(i,j,k) = drho2(i,j,k) + rho_int(i,j,k)*u_x(i,j,k) + u_cur(i,j,k)*rho_x(i,j,k)

        dv = 0.5d0*(LY5-LY1)/(rho_int(i,j,k)*c) + LZ4 + u_cur(i,j,k)*v_x(i,j,k)

        du = LY3 + LZ3 + p_x(i,j,k)/rho_int(i,j,k) + u_cur(i,j,k)*u_x(i,j,k)

        dw = 0.5d0*(LZ5-LZ1)/(rho_int(i,j,k)*c) + LY4 + u_cur(i,j,k)*w_x(i,j,k)

        dp = (LY5 + LY1)*0.5d0 + (LZ5 + LZ1)*0.5d0
        dp = dp + gama*p_int(i,j,k)*u_x(i,j,k) + u_cur(i,j,k)*p_x(i,j,k)

        drhou2(i,j,k) = u_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*du

        drhov2(i,j,k) = v_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dv

        drhow2(i,j,k) = w_cur(i,j,k)*drho2(i,j,k) + rho_int(i,j,k)*dw

        v_sqr = (u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2) * 0.5d0
        intr  = (u_cur(i,j,k)*du + v_cur(i,j,k)*dv + w_cur(i,j,k)*dw)*rho_int(i,j,k)
        de2(i,j,k) = v_sqr*drho2(i,j,k) + intr + dp/(gama-1.0d0)


        end do



        END SUBROUTINE
        !______________________________________________________________________


