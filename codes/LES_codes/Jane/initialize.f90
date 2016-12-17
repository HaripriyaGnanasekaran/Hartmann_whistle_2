        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 13-10-2009        Last Modified on 10-05-2010

        !_____________________________________________________________________


        
        !_____________________________________________________________________
        SUBROUTINE INITIALIZATION

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NY1,NY2
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dy1,dz1,jet_r
        USE FLOW_PARAMETERS,     ONLY : rho_atm,p_atm,v_jet,p_in,rho_in
        USE FLOW_PARAMETERS,     ONLY : c,gama,pi,p_c,rho_c,R
        USE continuity,          ONLY : rho_cur
        USE x_momentum,          ONLY : u_cur,rhou_cur
        USE y_momentum,          ONLY : v_cur,rhov_cur
        USE z_momentum,          ONLY : w_cur,rhow_cur
        USE pressure,            ONLY : p_cur
        USE energy,              ONLY : e_cur,T
        !USE vorticity,           ONLY : omega

        IMPLICIT NONE
       
        REAL (KIND = 8) :: x,y,z,rad
        REAL (KIND = 8) :: ke,M,d
        REAL (KIND = 8) :: rho_fac


        INTEGER :: i,j,k

        c = sqrt(1.4d0*101325.0d0/1.25d0)

        M = v_jet/c

        !x = 0.0d0
        !do i = NXL,NXR
        !y = 0.0d0
        !do j = NYB,NYT
        !z = 0.0d0
        !do k = NZB,NZF

        !rho_cur(i,j,k)  = rho_atm        
        !u_cur(i,j,k)    = 0.0D0
        !rhou_cur(i,j,k) = rho_cur(i,j,k)*u_cur(i,j,k)
        !v_cur(i,j,k)    = 0.0D0
        !rhov_cur(i,j,k) = rho_cur(i,j,k)*v_cur(i,j,k)
        !w_cur(i,j,k)    = 0.0d0
        !rhow_cur(i,j,k) = rho_cur(i,j,k)*w_cur(i,j,k)
        !p_cur(i,j,k)    = p_atm

        !ke              = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        !ke              = 0.5d0*rho_cur(i,j,k)*ke
        !e_cur(i,j,k)    = p_cur(i,j,k)/(gama-1.0d0) + ke
        !omega(i,j,k)    = 0.0d0

        !z = z + dz1
        !end do
        !y = y + dy1
        !end do  
        !x = x + dx1
        !end do
        
        !write(*,*) 'Main Initialization'
        d = (NY2-NY1)*dy1*0.5d0

        do k = NZB,NZF
        do j = NYB,NYT
        
        !y = (j - (NY1+NY2)/2)*dy1
        CALL Z_COORD(z,k)
        CALL Y_COORD(y,j)
        rad = sqrt(z*z+y*y)

        do i = NXL,NXR
        
        u_cur(i,j,k)      = v_jet*(0.5d0 + 0.5d0 * tanh((jet_r-rad)/(jet_r/10.0d0)) )

        if ( rad<= jet_r) then

        rho_cur(i,j,k)  = rho_in
        p_cur(i,j,k)    = p_in
        T(i,j,k)  	= p_cur(i,j,k)/(rho_cur(i,j,k) * R)


        else

        rho_cur(i,j,k)  = rho_c
        p_cur(i,j,k)    = p_c
        T(i,j,k)  	= p_cur(i,j,k)/(rho_cur(i,j,k) * R)

        end if

        rhou_cur(i,j,k) = rho_cur(i,j,k)*u_cur(i,j,k)        
        v_cur(i,j,k)    = 0.0D0
        rhov_cur(i,j,k) = rho_cur(i,j,k)*v_cur(i,j,k)
        w_cur(i,j,k)    = 0.0d0
        rhow_cur(i,j,k) = rho_cur(i,j,k)*w_cur(i,j,k)

        ke           = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke           = 0.5d0 * rho_cur(i,j,k)*ke
        e_cur(i,j,k) = p_cur(i,j,k)/(gama-1.0d0) + ke
        !omega(i,j,k) = 0.0d0

        end do
        end do        
        end do

        !write(*,*) 'Main initialization done'
        END SUBROUTINE
        !_____________________________________________________________________




        !_____________________________________________________________________
        SUBROUTINE DENSITY_INITIALIZATION

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF,NY1,NY2
        USE GRID_SIZE_TIME_STEP, ONLY : dz1,dz
        USE continuity,      ONLY : drho1,drho2,rho_x,rho_y,rho_z
        USE continuity,      ONLY : !rho_mean !,rho_ms
        USE continuity,      ONLY : rho_tot ,rho_ts
        USE continuity,      ONLY : rho_if
        USE FLOW_PARAMETERS, ONLY : rho_atm,pi

        IMPLICIT NONE
       
        INTEGER :: i,j,k
        REAL (KIND=8) :: zf,zb,z

        zf = NZB*dz1
        zb = NZF*dz1
      
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        drho1(i,j,k)    = 0.0d0
        drho2(i,j,k)    = 0.0d0
        rho_x(i,j,k)    = 0.0d0
        rho_y(i,j,k)    = 0.0d0
        rho_z(i,j,k)    = 0.0d0
        !rho_mean(i,j,k) = 0.0d0
        !rho_ms(i,j,k)   = 0.0d0
        rho_tot(i,j,k)  = 0.0d0
        rho_ts(i,j,k)   = 0.0d0


        end do
        end do
        end do

    
        END SUBROUTINE
        !_____________________________________________________________________



        !_____________________________________________________________________
        SUBROUTINE U_VELOCITY_INITIALIZATION

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE x_momentum,      ONLY : drhou1, drhou2,rhou_x!,rhou_y,rhou_z
        USE x_momentum,      ONLY : !u_mean !,u_ms
        USE x_momentum,      ONLY : u_tot ,u_ts,rhou_tot

        IMPLICIT NONE
       
        INTEGER :: i,j,k

      
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        drhou1(i,j,k) = 0.0d0
        drhou2(i,j,k) = 0.0d0
        rhou_x(i,j,k) = 0.0d0
        !rhou_y(i,j,k) = 0.0d0
        !rhou_z(i,j,k) = 0.0d0
        !u_mean(i,j,k) = 0.0d0
        !u_ms(i,j,k)   = 0.0d0
        u_tot(i,j,k)  = 0.0d0
        u_ts(i,j,k)   = 0.0d0
        rhou_tot(i,j,k)  = 0.0d0

        end do
        end do
        end do

        END SUBROUTINE
        !_____________________________________________________________________


        !_____________________________________________________________________
        SUBROUTINE V_VELOCITY_INITIALIZATION

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE y_momentum,  ONLY : drhov1, drhov2,rhov_y!,rhov_x,rhov_z
        USE y_momentum,  ONLY : !v_mean !,v_ms
        USE y_momentum,  ONLY : v_tot ,v_ts,rhov_tot

        IMPLICIT NONE
       
        INTEGER :: i,j,k


        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        drhov1(i,j,k) = 0.0d0
        drhov2(i,j,k) = 0.0d0
        !rhov_x(i,j,k) = 0.0d0
        rhov_y(i,j,k) = 0.0d0
        !rhov_z(i,j,k) = 0.0d0
        !v_mean(i,j,k) = 0.0d0
        !v_ms(i,j,k)   = 0.0d0
        v_tot(i,j,k)  = 0.0d0
        v_ts(i,j,k)   = 0.0d0
        rhov_tot(i,j,k)  = 0.0d0

        end do
        end do
        end do


        END SUBROUTINE
        !_____________________________________________________________________



        !_____________________________________________________________________
        SUBROUTINE W_VELOCITY_INITIALIZATION

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE z_momentum,  ONLY : drhow1, drhow2,rhow_z!,rhow_y,rhow_x
        USE z_momentum,  ONLY : !w_mean !,w_ms
        USE z_momentum,  ONLY : w_tot ,w_ts,rhow_tot

        IMPLICIT NONE
       
        INTEGER :: i,j,k


        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        drhow1(i,j,k) = 0.0d0
        drhow2(i,j,k) = 0.0d0
        !rhow_x(i,j,k) = 0.0d0
        !rhow_y(i,j,k) = 0.0d0
        rhow_z(i,j,k) = 0.0d0
        !w_mean(i,j,k) = 0.0d0
        !w_ms(i,j,k)   = 0.0d0
        w_tot(i,j,k)  = 0.0d0
        w_ts(i,j,k)   = 0.0d0
        rhow_tot(i,j,k)  = 0.0d0

        end do
        end do
        end do


        END SUBROUTINE
        !_____________________________________________________________________



        !_____________________________________________________________________
        SUBROUTINE PRESSURE_INITIALIZATION

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE pressure,  ONLY : p_x,p_y,p_z !,p_mean,p_ms
        USE pressure,  ONLY : p_tot,p_ts

        IMPLICIT NONE
       
        INTEGER :: i,j,k

      
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        p_x(i,j,k)    = 0.0d0
        p_y(i,j,k)    = 0.0d0
        p_z(i,j,k)    = 0.0d0
        !p_mean(i,j,k) = 0.0d0
        !p_ms(i,j,k)   = 0.0d0
        p_tot(i,j,k)  = 0.0d0
        p_ts(i,j,k)   = 0.0d0

        end do
        end do
        end do

        END SUBROUTINE
        !_____________________________________________________________________


        !_____________________________________________________________________
        SUBROUTINE ENERGY_INITIALIZATION

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE energy,          ONLY : de1,de2 !,e_x,e_y,e_z
        USE energy,          ONLY : e_tot !,e_mean !e_ms,e_ts
        USE energy,          ONLY : vke_xx,vke_yy,vke_zz

        IMPLICIT NONE
       
        INTEGER :: i,j,k

      
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        de1(i,j,k)     = 0.0d0
        de2(i,j,k)     = 0.0d0
        !e_x(i,j,k)     = 0.0d0
        !e_y(i,j,k)     = 0.0d0
        !e_z(i,j,k)     = 0.0d0
        !e_mean(i,j,k)  = 0.0d0
        !e_ms(i,j,k)    = 0.0d0
        vke_xx(i,j,k)   = 0.0d0
        vke_yy(i,j,k)   = 0.0d0
        vke_zz(i,j,k)   = 0.0d0
        e_tot(i,j,k)   = 0.0d0
        !e_ts(i,j,k)    = 0.0d0

        end do
        end do
        end do

        END SUBROUTINE
        !_____________________________________________________________________




        !_____________________________________________________________________
        SUBROUTINE VORTICITY_INITIALIZATION

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE vorticity,       ONLY : omegax,omegay,omegaz
        !USE vorticity,       ONLY : omega_tot,omega_ts

        IMPLICIT NONE
       
        INTEGER :: i,j,k

      
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        omegax(i,j,k)     = 0.0d0
        omegay(i,j,k)     = 0.0d0
        omegaz(i,j,k)     = 0.0d0
        !omega_tot(i,j,k)  = 0.0d0
        !omega_ts(i,j,k)   = 0.0d0

        end do
        end do
        end do

        END SUBROUTINE
        !_____________________________________________________________________



        !_____________________________________________________________________
        SUBROUTINE CROSS_MOMENTUM_INITIALIZATION

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE cross_momentum , ONLY : uv_tot,uw_tot,vw_tot,turb_ke_tot

        IMPLICIT NONE
       
        INTEGER :: i,j,k

      
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        uv_tot(i,j,k) = 0.0d0
        uw_tot(i,j,k) = 0.0d0
        vw_tot(i,j,k) = 0.0d0
        turb_ke_tot(i,j,k) = 0.0d0
        end do
        end do
        end do

        END SUBROUTINE
        !_____________________________________________________________________






        !_____________________________________________________________________
        SUBROUTINE LIGHTHILL_INITIALIZATION

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE lighthill_source

        IMPLICIT NONE
       
        INTEGER :: i,j,k

      
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        lh_s(i,j,k)      = 0.0d0
        lh_s_tot(i,j,k)  = 0.0d0
        lh_s_ts(i,j,k)   = 0.0d0

        lh_s_f(i,j,k)    = 0.0d0
        lh_s_b(i,j,k)    = 0.0d0


        end do
        end do
        end do

        END SUBROUTINE
        !_____________________________________________________________________





        !_____________________________________________________________________
        SUBROUTINE GRID_INITIALIZATION

        USE GRID_DIMENSIONS
        USE GRID_SIZE_TIME_STEP


        IMPLICIT NONE
       
        INTEGER :: i,j,k

        !do k = NZB,NZF
        !do j = NYB,NYT
        do i = NXL,NX_st
        dxi(i) = 1.0d0/dx1
        end do
        
        do i=NX_st,NXO-1
        dxi(i) = dxi(i-1)/axs
        end do
        
        do i = NXO,NXR-1
        dxi(i) = dxi(i-1)/axo
        end do
        !end do
        !end do


        !do k = NZB,NZF
        !do i = NXL,NXR
        dyi(NYB) = 1.0d0/(dy1*ayo**(NYSB-NYB))
        do j = NYB+1,NYSB-1
        dyi(j) = dyi(j-1)*ayo
        end do
        
        do j = NYSB,NYST-1
        dyi(j) = 1.0d0/dy1
        end do

        do j = NYST,NYT-1
        dyi(j) = dyi(j-1)/ayo
        end do
        !end do
        !end do
        

        !do k = NZB,NZF
        !do i = NXL,NXR
        dzi(NZB) = 1.0d0/(dz1*azo**(NZSB-NZB))
        do k = NZB+1,NZSF-1
        dzi(k) = dzi(k-1)*azo
        end do
        
        do k = NZSB,NZSF-1
        dzi(k) = 1.0d0/dz1
        end do

        do k = NZSF,NZF-1
        dzi(k) = dzi(k-1)/azo
        end do
        !end do
        !end do



        END SUBROUTINE
        !_____________________________________________________________________
	SUBROUTINE METRICS_X

	USE GRID_DIMENSIONS
	USE GRID_SIZE_TIME_STEP
        USE METRICS
        	
	REAL(KIND=8) :: d_eta,d_eta_inv,temp,a,adi
	REAL(KIND=8),DIMENSION(NXL:NXR), SAVE :: funct,x,dx_deta_f,dx_deta_b
	a  =  0.21132486540518708
        adi =  0.78867513459481287
	d_eta = dy1
	d_eta_inv = 1.0d0/d_eta
	
	DO i=NXL,NXR
		CALL X_COORD(temp,i)
		x(i) = temp
	END DO
	!*****************************************
	!COMPACT HIXON FORWARD SCHEME TO DETERMINE DX/DE 
	!
	DO i= NXL,NXR-1
		funct(i) = (x(i+1)-x(i))*d_eta_inv
	END DO
	!BOUNDARY FORMULATION
	dx_deta_f(NXR) = -4.0d0*(X(NXR-1)) + X(NXR-2) + 3.0D0*(X(NXR))
	dx_deta_f(NXR) = 0.5d0*dx_deta_f(NXR)*d_eta_inv

	DO i = NXR-1,NXL,-1
		dx_deta_f(i) = (funct(i) - a*(dx_deta_f(i+1)))*adi
	END DO

      
	!******************************************
	!*****************************************
	!COMPACT HIXON BACKWARD SCHEME TO DETERMINE DX/DE 
	!
        
	!BOUNDARY FORMULATION
	dx_deta_b(NXL) = 4.0d0*X(NXL+1)-X(NXL+2)-3.0D0*X(NXL)
        dx_deta_b(NXL) = 0.5d0*dx_deta_b(NXL)*d_eta_inv
        
	DO i = NXL,NXR-1
		dx_deta_b(i+1) = (funct(i) - a*(dx_deta_b(i)))*adi
	END DO
	!******************************************	
	DO i = NXL,NXR
		dx_deta(i) = 0.5d0*(dx_deta_f(i)+dx_deta_b(i))
	END DO
	
        END SUBROUTINE
        !____________________________________________________________________
        SUBROUTINE METRICS_Y
        USE GRID_DIMENSIONS
        USE GRID_SIZE_TIME_STEP
        USE METRICS

        INTEGER :: i
        REAL(KIND=8) :: d_eta,d_eta_inv,temp,a,adi
        REAL(KIND=8),DIMENSION(NYB:NYT), SAVE ::funct,y,dy_deta_f,dy_deta_b
        a  =  0.21132486540518708
        adi =  0.78867513459481287
        d_eta = dy1
        d_eta_inv = 1.0d0/d_eta
        
        DO i=NYB,NYT
                CALL Y_COORD(temp,i)
                y(i) = temp
        END DO
        !*****************************************
        !COMPACT HIXON FORWARD SCHEME TO DETERMINE DX/DE 
        !
        DO i= NYB,NYT-1
                funct(i) = (y(i+1)-y(i))*d_eta_inv
        END DO
        !BOUNDARY FORMULATION
        dy_deta_f(NYT) = -4.0d0*(y(NYT-1)) + y(NYT-2) + 3.0D0*(y(NYT))
        dy_deta_f(NYT) = 0.5d0*dy_deta_f(NYT)*d_eta_inv

        DO i = NYT-1,NYB,-1
                dy_deta_f(i) = (funct(i) - a*(dy_deta_f(i+1)))*adi
        END DO


        !******************************************
        !*****************************************
        !COMPACT HIXON BACKWARD SCHEME TO DETERMINE DX/DE 
        !
        
        !BOUNDARY FORMULATION
        dy_deta_b(NYB) = 4.0d0*y(NYB+1)-y(NYB+2)-3.0D0*y(NYB)
        dy_deta_b(NYB) = 0.5d0*dy_deta_b(NYB)*d_eta_inv

        DO i = NYB,NYT-1
                dy_deta_b(i+1) = (funct(i) - a*(dy_deta_b(i)))*adi
        END DO
        !******************************************     
        DO i = NYB,NYT
                dy_deta(i) = 0.5d0*(dy_deta_f(i)+dy_deta_b(i))
        END DO
        END SUBROUTINE
        !_________________________________________________________________
        !____________________________________________________________________
        SUBROUTINE METRICS_Z
        USE GRID_DIMENSIONS
        USE GRID_SIZE_TIME_STEP
        USE METRICS

        INTEGER :: i
        REAL(KIND=8) :: d_eta,d_eta_inv,temp,a,adi
        REAL(KIND=8),DIMENSION(NZB:NZF), SAVE::funct,z,dz_deta_f,dz_deta_b
        a  =  0.21132486540518708
        adi =  0.78867513459481287
        d_eta = dy1
        d_eta_inv = 1.0d0/d_eta

        DO i=NZB,NZF
                CALL Z_COORD(temp,i)
                z(i) = temp
        END DO
        !*****************************************
        !COMPACT HIXON FORWARD SCHEME TO DETERMINE DX/DE 
        !
        DO i= NZB,NZF-1
                funct(i) = (z(i+1)-z(i))*d_eta_inv
        END DO
        !BOUNDARY FORMULATION
        dz_deta_f(NZF) = -4.0d0*(z(NZF-1)) + z(NZF-2) + 3.0D0*(z(NZF))
        dz_deta_f(NZF) = 0.5d0*dz_deta_f(NZF)*d_eta_inv

        DO i = NZF-1,NZB,-1
                dz_deta_f(i) = (funct(i) - a*(dz_deta_f(i+1)))*adi
        END DO


        !******************************************
        !*****************************************
        !COMPACT HIXON BACKWARD SCHEME TO DETERMINE DX/DE 
        !

        !BOUNDARY FORMULATION
        dz_deta_b(NZB) = 4.0d0*z(NZB+1)-z(NZB+2)-3.0D0*z(NZB)
        dz_deta_b(NZB) = 0.5d0*dz_deta_b(NZB)*d_eta_inv

        DO i = NZB,NZF-1
                dz_deta_b(i+1) = (funct(i) - a*(dz_deta_b(i)))*adi
        END DO
        !******************************************     
        DO i = NZB,NZF
               dz_deta(i) = 0.5d0*(dz_deta_f(i)+dz_deta_b(i))
        END DO

        END SUBROUTINE
        !_________________________________________________________________



        !_____________________________________________________________________
        SUBROUTINE INITIALIZE_RANDOM_NUMBERS_FREUND
        USE turbulence_freund
        USE FLOW_PARAMETERS, ONLY : v_jet,pi
        IMPLICIT NONE 
        INTEGER :: i,j,k,l,m


	dA	= 0.0001d0
	d_alpha = 0.00085d0

        CALL random_number(A)	
        CALL random_number(alpha)
        CALL random_number(phi)
        CALL random_number(psi)

	do l = 1,NTF
	do m = 0,NTA-1
	A(l,m) = 0.01d0 + 0.06d0*A(l,m)
	alpha(l,m) = 0.1d0 + 0.6d0*alpha(l,m)
        phi(l,m) = phi(l,m) * 2.0d0 * pi
        psi(l,m) = psi(l,m) * 2.0d0 * pi
	end do
	end do


 
        END SUBROUTINE
        !_____________________________________________________________________



        !_____________________________________________________________________
        SUBROUTINE INFLOW_FREUND(ts)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF,NY1,NY2
        USE GRID_SIZE_TIME_STEP, ONLY : dz1,dz,dt
        USE GRID_SIZE_TIME_STEP, ONLY : jet_d,jet_r
        USE FLOW_PARAMETERS, ONLY : pi,v_jet
        USE continuity,      ONLY : rho_if
        USE x_momentum,      ONLY : u_tar
        USE turbulence_freund

        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ts       
        INTEGER :: i,j,k,l,m
        REAL (KIND=8) :: z,t,d,x,y,r,r_nd,r_exp
        REAL (KIND=8) :: v_theta,theta,v_r
        REAL (KIND=8) :: R_turb,char_freq
	REAL (KIND=8) :: b,dif


        t  = dt * ts

        char_freq =2.0d0 * pi *  v_jet/jet_d !* 20.0d0

	do l = 1,NTF
	do m = 0,NTA-1
	CALL random_number(dif)
	if (dif < 1.0d0/3.0d0) then
	dif = -1.0d0
	else if ((dif >= 1.0d0/3.0d0) .and. (dif < 2.0d0/3.0d0 ) ) then
	dif = 0.0d0
	else
	dif = 1.0d0
	end if

	if( ((A(l,m)+dif*dA) >= 0.01d0) .and. ((A(l,m)+dif*dA) <= 0.07d0)  ) then
	A(l,m) = A(l,m)+dif*dA
	end if

	if( ((alpha(l,m)+dif*d_alpha) >= 0.1d0) .and. ((alpha(l,m)+dif*d_alpha) <= 0.7d0)  ) then
	alpha(l,m) = alpha(l,m)+dif*d_alpha
	end if

	end do
	end do

        do k = NZB,NZF
        do j = NYB,NYT

        u_tar(j,k) = 0.0d0

        CALL y_coord(y,j)
        CALL z_coord(z,k)

        r = sqrt(y*y+z*z)
        if (r .ne. 0.0d0) then
        theta = atan2(z,y)
        end if
	b= 10.0d0
        if (r<=jet_r) then

        do l = 0,NTF
        do m = 0,NTA-1
        b = b + A(l,m)*cos(char_freq*alpha(l,m)*t+phi(l,m))*cos(m*theta+psi(l,m))        
        end do
	end do


        u_tar(j,k) = v_jet * (0.5d0 + 0.5d0*tanh((jet_r-r)/(jet_r/b)) )  
	
	end if

        end do
        end do
        !end if



        END SUBROUTINE
        !_____________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE WRITE_INFLOW_FREUND(l)
        USE turbulence_freund    
      

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 40) :: fileoutput

        INTEGER :: i,j,k


        fileoutput = 'turbulence_freund_000000.dat'

        WRITE(fileoutput,'(a18,I6.6,a)') fileoutput,l,'.dat'

        OPEN(UNIT = 144, FILE = fileoutput , ACTION ='WRITE')

	do i = 1,NTF
	do j = 0,NTA-1
	write(144,*) A(i,j), alpha(i,j),phi(i,j),psi(i,j)
        end do
        end do


        CLOSE ( UNIT = 144)


        END SUBROUTINE
        !_____________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE READ_INFLOW_FREUND(l)
        USE turbulence_freund    
      

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 40) :: fileoutput

        INTEGER :: i,j,k


        fileoutput = 'turbulence_freund_000000.dat'

        WRITE(fileoutput,'(a18,I6.6,a)') fileoutput,l,'.dat'

        OPEN(UNIT = 144, FILE = fileoutput , ACTION ='READ')

	do i = 1,NTF
	do j = 0,NTA-1
	read(144,*) A(i,j), alpha(i,j),phi(i,j),psi(i,j)
        end do
        end do


        CLOSE ( UNIT = 144)


        END SUBROUTINE
        !_____________________________________________________________________



