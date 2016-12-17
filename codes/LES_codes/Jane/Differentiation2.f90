        !______________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 11-07-2009        Last Modified on 16-05-2010

        !______________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE EVALUATION2(ch)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_SIZE_TIME_STEP, ONLY : dt

        USE continuity, ONLY : rho_int,drho2,rho_x,rho_y,rho_z

        USE x_momentum, ONLY : rhou_int,drhou2,rhou_x!,rhou_y,rhou_z
        USE x_momentum, ONLY : vrhou_xx,vrhou_yy,vrhou_zz
        USE x_momentum, ONLY : u_cur,u_x,u_y,u_z

        USE y_momentum, ONLY : rhov_int,drhov2,rhov_y!,rhov_x,rhov_z
        USE y_momentum, ONLY : vrhov_xx,vrhov_yy,vrhov_zz
        USE y_momentum, ONLY : v_cur,v_x,v_y,v_z

        USE z_momentum, ONLY : rhow_int,drhow2,rhow_z!,rhow_y,rhow_x
        USE z_momentum, ONLY : vrhow_xx,vrhow_yy,vrhow_zz
        USE z_momentum, ONLY : w_cur,w_x,w_y,w_z

        USE pressure,   ONLY : p_int,p_x,p_y,p_z

        USE energy,     ONLY : e_int,de2,vke_xx,vke_yy,vke_zz
        USE energy,     ONLY : efl_x,efl_y,efl_z,T,T_x,T_y,T_z

        USE x_momentum, ONLY : rhouu_x,rhouv_y,rhouw_z
        USE y_momentum, ONLY : rhovu_x,rhovv_y,rhovw_z
        USE z_momentum, ONLY : rhowu_x,rhowv_y,rhoww_z

        USE FLOW_PARAMETERS, ONLY : gama,c,mu,kT,Cp,Pr

        IMPLICIT NONE
        
        INTEGER, DIMENSION(2),              INTENT(IN) :: ch

        INTEGER :: i,j,k
	REAL (KIND=8) :: m1,m2,m3
        m1 = mu
	m2 = mu*4.0d0/3.0d0
	m3 = mu*2.0d0/3.0d0

        !c = sqrt(gama*p_atm/rho_atm)

        !write(*,*) 'Differentiation'

                !!CONTINUITY EQUATION
                !rhou_x = rhou_int !CALL COPY(rhou_int,rhou_x)
                !rhov_y = rhov_int !CALL COPY(rhov_int,rhov_y)
                !rhow_z = rhow_int !CALL COPY(rhow_int,rhow_z)
                !!CALL COPY(rho_int,rho_x)
                !!CALL COPY(rho_int,rho_y)
                !!CALL COPY(rho_int,rho_z)

                !!MOMENTUM EQUATION VISCOUS TERMS
                !u_x = u_cur !CALL COPY(u_cur,u_x)
                !u_y = u_cur !CALL COPY(u_cur,u_y)
                !u_z = u_cur !CALL COPY(u_cur,u_z)
                !v_x = v_cur !CALL COPY(v_cur,v_x)
                !v_y = v_cur !CALL COPY(v_cur,v_y)
                !v_z = v_cur !CALL COPY(v_cur,v_z)
                !w_x = w_cur !CALL COPY(w_cur,w_x)
                !w_y = w_cur !CALL COPY(w_cur,w_y)
                !w_z = w_cur !CALL COPY(w_cur,w_z)

        !X_MOMENTUM CONVECTIVE TERMS
        CALL PRODUCT_OF_TWO_TERMS(rhou_int,u_cur,rhouu_x)
        !CALL DIVIDE(rho_int,rhouu_x)
        CALL PRODUCT_OF_TWO_TERMS(rhou_int,v_cur,rhouv_y)
        !CALL DIVIDE(rho_int,rhouv_y)
        CALL PRODUCT_OF_TWO_TERMS(rhou_int,w_cur,rhouw_z)
        !CALL DIVIDE(rho_int,rhouw_z)
        CALL COPY(p_int,p_x)

        !Y_MOMENTUM CONVECTIVE TERMS
        CALL COPY(rhouv_y,rhovu_x)
        CALL PRODUCT_OF_TWO_TERMS(rhov_int,v_cur,rhovv_y)
        !CALL DIVIDE(rho_int,rhovv_y)
        CALL PRODUCT_OF_TWO_TERMS(rhov_int,w_cur,rhovw_z)
        !CALL DIVIDE(rho_int,rhovw_z)
        CALL COPY(p_int,p_y)

        !Z_MOMENTUM CONVECTIVE TERMS
        CALL COPY(rhouw_z,rhowu_x)
        CALL COPY(rhovw_z,rhowv_y)
        CALL PRODUCT_OF_TWO_TERMS(rhow_int,w_cur,rhoww_z)
        !CALL DIVIDE(rho_int,rhoww_z)
        CALL COPY(p_int,p_z)

        !ENERGY EQUATION CONVECTIVE TERMS
        CALL SUM_OF_TWO_TERMS(e_int,p_int,efl_x)
        CALL COPY(efl_x,efl_y)
        CALL COPY(efl_x,efl_z)
        CALL MULTIPLY(u_cur,efl_x)
        CALL MULTIPLY(v_cur,efl_y)
        CALL MULTIPLY(w_cur,efl_z)

        !write(*,*) 'Differentiation copy 1'



        if (ch(1) == 0) then
        !CONTINUITY EQUATION
        CALL TURKEL_FORWARD_X(rhou_int,rhou_x)
        CALL TURKEL_FORWARD_X(rho_int,rho_x)
        !MOMENTUM EQUATION VISCOUS TERMS
        CALL TURKEL_FORWARD_X(u_cur,u_x)
        CALL TURKEL_FORWARD_X(v_cur,v_x)
        CALL TURKEL_FORWARD_X(w_cur,w_x)
        !X_MOMENTUM CONVECTIVE TERM
        CALL HIXON_FORWARD_X(rhouu_x)
        CALL HIXON_FORWARD_X(p_x)
        !Y_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_FORWARD_X(rhovu_x)
        !Z_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_FORWARD_X(rhowu_x)
        !ENERGY EQUATION CONVECTIVE TERMS
        CALL HIXON_FORWARD_X(efl_x)
        !ENERGY EQUATION CONDUCTION TERM
        CALL TURKEL_FORWARD_X(T,T_x)


        !CONTINUITY EQUATION
        CALL TURKEL_FORWARD_Y(rhov_int,rhov_y)
        CALL TURKEL_FORWARD_Y(rho_int,rho_y)
        !MOMENTUM EQUATION VISCOUS TERMS
        CALL TURKEL_FORWARD_Y(u_cur,u_y)
        CALL TURKEL_FORWARD_Y(v_cur,v_y)
        CALL TURKEL_FORWARD_Y(w_cur,w_y)
        !X_MOMENTUM CONVECTIVE TERM
        CALL HIXON_FORWARD_Y(rhouv_y)
        !Y_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_FORWARD_Y(rhovv_y)
        CALL HIXON_FORWARD_Y(p_y)
        !Z_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_FORWARD_Y(rhowv_y)
        !ENERGY EQUATION CONVECTIVE TERMS
        CALL HIXON_FORWARD_Y(efl_y)
        !ENERGY EQUATION CONDUCTION TERM
        CALL TURKEL_FORWARD_Y(T,T_y)


        !CONTINUITY EQUATION
        CALL TURKEL_FORWARD_Z(rhow_int,rhow_z)
        CALL TURKEL_FORWARD_Z(rho_int,rho_z)
        !MOMENTUM EQUATION VISCOUS TERMS
        CALL TURKEL_FORWARD_Z(u_cur,u_z)
        CALL TURKEL_FORWARD_Z(v_cur,v_z)
        CALL TURKEL_FORWARD_Z(w_cur,w_z)
        !X_MOMENTUM CONVECTIVE TERM
        CALL HIXON_FORWARD_Z(rhouw_z)
        !Y_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_FORWARD_Z(rhovw_z)
        !Z_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_FORWARD_Z(rhoww_z)
        CALL HIXON_FORWARD_Z(p_z)
        !ENERGY EQUATION CONVECTIVE TERMS
        CALL HIXON_FORWARD_Z(efl_z)
        !ENERGY EQUATION CONDUCTION TERM
        CALL TURKEL_FORWARD_Z(T,T_z)


        else


        !CONTINUITY EQUATION
        CALL TURKEL_BACKDIF_X(rhou_int,rhou_x)
        CALL TURKEL_BACKDIF_X(rho_int,rho_x)
        !MOMENTUM EQUATION VISCOUS TERMS
        CALL TURKEL_BACKDIF_X(u_cur,u_x)
        CALL TURKEL_BACKDIF_X(v_cur,v_x)
        CALL TURKEL_BACKDIF_X(w_cur,w_x)
        !X_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_BACKDIF_X(rhouu_x)
        CALL HIXON_BACKDIF_X(p_x)
        !Y_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_BACKDIF_X(rhovu_x)
        !Z_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_BACKDIF_X(rhowu_x)
        !ENERGY EQUATION CONVECTIVE TERMS
        CALL HIXON_BACKDIF_X(efl_x)
        !ENERGY EQUATION CONDUCTION TERM
        CALL TURKEL_BACKDIF_X(T,T_x)


        !CONTINUITY EQUATION
        CALL TURKEL_BACKDIF_Y(rhov_int,rhov_y)
        CALL TURKEL_BACKDIF_Y(rho_int,rho_y)
        !MOMENTUM EQUATION VISCOUS TERMS
        CALL TURKEL_BACKDIF_Y(u_cur,u_y)
        CALL TURKEL_BACKDIF_Y(v_cur,v_y)
        CALL TURKEL_BACKDIF_Y(w_cur,w_y)
        !X_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_BACKDIF_Y(rhouv_y)
        !Y_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_BACKDIF_Y(rhovv_y)
        CALL HIXON_BACKDIF_Y(p_y)
        !Z_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_BACKDIF_Y(rhowv_y)
        !ENERGY EQUATION CONVECTIVE TERMS
        CALL HIXON_BACKDIF_Y(efl_y)
        !ENERGY EQUATION CONDUCTION TERM
        CALL TURKEL_BACKDIF_Y(T,T_y)


        !CONTINUITY EQUATION
        CALL TURKEL_BACKDIF_Z(rhow_int,rhow_z)
        CALL TURKEL_BACKDIF_Z(rho_int,rho_z)
        !MOMENTUM EQUATION VISCOUS TERMS
        CALL TURKEL_BACKDIF_Z(u_cur,u_z)
        CALL TURKEL_BACKDIF_Z(v_cur,v_z)
        CALL TURKEL_BACKDIF_Z(w_cur,w_z)   
        !X_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_BACKDIF_Z(rhouw_z)
        !Y_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_BACKDIF_Z(rhovw_z)
        !Z_MOMENTUM CONVECTIVE TERMS
        CALL HIXON_BACKDIF_Z(rhoww_z)
        CALL HIXON_BACKDIF_Z(p_z)     
        !ENERGY EQUATION CONVECTIVE TERMS
        CALL HIXON_BACKDIF_Z(efl_z)
        !ENERGY EQUATION CONDUCTION TERM
        CALL TURKEL_BACKDIF_Z(T,T_z)


        end if


        !ALLOCATE(vrhou_xx(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(vrhov_xx(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(vrhow_xx(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(vrhou_yy(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(vrhov_yy(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(vrhow_yy(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(vrhou_zz(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(vrhov_zz(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(vrhow_zz(NXL:NXR,NYB:NYT,NZB:NZF))

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !vrhou_xx(i,j,k) = 0.0d0
        !vrhov_xx(i,j,k) = 0.0d0
        !vrhow_xx(i,j,k) = 0.0d0
        !vrhou_yy(i,j,k) = 0.0d0
        !vrhov_yy(i,j,k) = 0.0d0
        !vrhow_yy(i,j,k) = 0.0d0
        !vrhou_zz(i,j,k) = 0.0d0
        !vrhov_zz(i,j,k) = 0.0d0
        !vrhow_zz(i,j,k) = 0.0d0
        !end do
        !end do
        !end do


        
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        CALL SUTHERLAND(T(i,j,k),mu)

        m1 = mu
	m2 = mu * 4.0d0/3.0d0
	m3 = mu * 2.0d0/3.0d0

        vrhou_xx(i,j,k) = m2*u_x(i,j,k) - m3*(v_y(i,j,k) + w_z(i,j,k))

        vrhou_yy(i,j,k) = mu*(u_y(i,j,k) + v_x(i,j,k))

        vrhou_zz(i,j,k) = mu*(u_z(i,j,k) + w_x(i,j,k))

        vrhov_xx(i,j,k) = vrhou_yy(i,j,k)

        vrhov_yy(i,j,k) = m2*v_y(i,j,k) - m3*(u_x(i,j,k) + w_z(i,j,k))

        vrhov_zz(i,j,k) = mu*(v_z(i,j,k) + w_y(i,j,k))

        vrhow_xx(i,j,k) = vrhou_zz(i,j,k)

        vrhow_yy(i,j,k) = vrhov_zz(i,j,k)

        vrhow_zz(i,j,k) = m2*w_z(i,j,k) - m3*(u_x(i,j,k) + v_y(i,j,k))

        kT = mu*Cp/Pr

        vke_xx(i,j,k) = vrhou_xx(i,j,k)*u_cur(i,j,k) + vrhou_yy(i,j,k)*v_cur(i,j,k)
        vke_xx(i,j,k) = vke_xx(i,j,k) + vrhou_zz(i,j,k)*w_cur(i,j,k) + kT*T_x(i,j,k)

        vke_yy(i,j,k) = vrhov_xx(i,j,k)*u_cur(i,j,k) + vrhov_yy(i,j,k)*v_cur(i,j,k)
        vke_yy(i,j,k) = vke_yy(i,j,k) + vrhov_zz(i,j,k)*w_cur(i,j,k) + kT*T_y(i,j,k)

        vke_zz(i,j,k) = vrhow_xx(i,j,k)*u_cur(i,j,k) + vrhow_yy(i,j,k)*v_cur(i,j,k)
        vke_zz(i,j,k) = vke_zz(i,j,k) + vrhow_zz(i,j,k)*w_cur(i,j,k) + kT*T_z(i,j,k)

        end do
        end do
        end do

        
        if  (ch(2) == 0) then
        CALL HIXON_FORWARD_X(vrhou_xx)
        CALL HIXON_FORWARD_X(vrhov_xx)
        CALL HIXON_FORWARD_X(vrhow_xx)
        CALL HIXON_FORWARD_X(vke_xx)

        CALL HIXON_FORWARD_Y(vrhou_yy)
        CALL HIXON_FORWARD_Y(vrhov_yy)
        CALL HIXON_FORWARD_Y(vrhow_yy)
        CALL HIXON_FORWARD_Y(vke_yy)

        CALL HIXON_FORWARD_Z(vrhou_zz)
        CALL HIXON_FORWARD_Z(vrhov_zz)
        CALL HIXON_FORWARD_Z(vrhow_zz)
        CALL HIXON_FORWARD_Z(vke_zz)

        else
        CALL HIXON_BACKDIF_X(vrhou_xx)
        CALL HIXON_BACKDIF_X(vrhov_xx)
        CALL HIXON_BACKDIF_X(vrhow_xx)
        CALL HIXON_BACKDIF_X(vke_xx)

        CALL HIXON_BACKDIF_Y(vrhou_yy)
        CALL HIXON_BACKDIF_Y(vrhov_yy)
        CALL HIXON_BACKDIF_Y(vrhow_yy)
        CALL HIXON_BACKDIF_Y(vke_yy)

        CALL HIXON_BACKDIF_Z(vrhou_zz)
        CALL HIXON_BACKDIF_Z(vrhov_zz)
        CALL HIXON_BACKDIF_Z(vrhow_zz)
        CALL HIXON_BACKDIF_Z(vke_zz)

        end if
       
        !write(*,*) 'Hixon Completed'
        
        do k = NZB+1,NZF-1
        do j = NYB+1,NYT-1
        do i = NXL+1,NXR-1

        drho2(i,j,k) = rhou_x(i,j,k) + rhov_y(i,j,k) + rhow_z(i,j,k)

        drhou2(i,j,k) = rhouu_x(i,j,k) + rhouv_y(i,j,k) + rhouw_z(i,j,k)
        drhou2(i,j,k) = drhou2(i,j,k) + p_x(i,j,k)
        drhou2(i,j,k) = drhou2(i,j,k) - (vrhou_xx(i,j,k)+vrhou_yy(i,j,k)+vrhou_zz(i,j,k))
      
        drhov2(i,j,k) = rhovu_x(i,j,k) + rhovv_y(i,j,k) + rhovw_z(i,j,k)
        drhov2(i,j,k) = drhov2(i,j,k) + p_y(i,j,k)
        drhov2(i,j,k) = drhov2(i,j,k) - (vrhov_xx(i,j,k)+vrhov_yy(i,j,k)+vrhov_zz(i,j,k))

        drhow2(i,j,k) = rhowu_x(i,j,k) + rhowv_y(i,j,k) + rhoww_z(i,j,k)
        drhow2(i,j,k) = drhow2(i,j,k) + p_z(i,j,k)
        drhow2(i,j,k) = drhow2(i,j,k) - (vrhow_xx(i,j,k)+vrhow_yy(i,j,k)+vrhow_zz(i,j,k))

        de2(i,j,k) = efl_x(i,j,k) + efl_y(i,j,k) + efl_z(i,j,k)
        de2(i,j,k) = de2(i,j,k) - (vke_xx(i,j,k)+vke_yy(i,j,k)+vke_zz(i,j,k))
  
        end do
        end do
        end do


        END SUBROUTINE
        !______________________________________________________________________




