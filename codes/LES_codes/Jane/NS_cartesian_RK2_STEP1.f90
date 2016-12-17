        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 11-07-2009        Last Modified on 30-05-2010

        !_____________________________________________________________________




        !_____________________________________________________________________
        SUBROUTINE RK2_STEP1
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NY1,NY2
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dy1,dz1,dt,jet_r

        USE FLOW_PARAMETERS

        USE continuity, ONLY : rho_cur,rho_int,drho1,rho_if
        USE x_momentum 
        USE y_momentum, ONLY : rhov_cur,rhov_int,drhov1,v_cur
        USE y_momentum, ONLY : v_if
        USE z_momentum, ONLY : rhow_cur,rhow_int,drhow1,w_cur
        USE z_momentum, ONLY : w_if
        USE pressure,   ONLY : p_cur,p_int
        USE energy,     ONLY : e_cur,e_int,de1
        IMPLICIT NONE

        INTEGER :: i,j,k

        REAL (KIND=8) :: ke,d,v_in,M,sig,sigi,u_in,rho_inlet
        REAL (KIND=8) :: x,y,z,rad,rho_fac 

        c = sqrt(gama*p_in/rho_in)
        M = v_jet/c
        d = (NY2-NY1)*dy1*0.5d0

        sig     = 0.005d0
        sigi    = 0.05d0


        !ALLOCATE(rho_int(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(rhou_int(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(rhov_int(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(rhow_int(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(p_int(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(e_int(NXL:NXR,NYB:NYT,NZB:NZF))

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF

        !rho_int(i,j,k)  = 0.0d0
        !rhou_int(i,j,k) = 0.0d0
        !rhov_int(i,j,k) = 0.0d0
        !rhow_int(i,j,k) = 0.0d0
        !p_int(i,j,k)    = 0.0d0
        !e_int(i,j,k)    = 0.0d0

        !end do
        !end do
        !end do


        !write(*,*) 'After Allocation STEP1' 
        !write(*,*) p_int(NXL,NYB,NZF),e_int(NXL,NYB,NZF)
        !write(*,*) p_cur(NXL,NYB,NZF),e_cur(NXL,NYB,NZF)

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        rho_int(i,j,k)  = rho_cur(i,j,k) - drho1(i,j,k)*dt

        rhou_int(i,j,k) = rhou_cur(i,j,k) - drhou1(i,j,k)*dt

        rhov_int(i,j,k) = rhov_cur(i,j,k) - drhov1(i,j,k)*dt

        rhow_int(i,j,k) = rhow_cur(i,j,k) - drhow1(i,j,k)*dt

        e_int(i,j,k)    = e_cur(i,j,k) - de1(i,j,k)*dt 

        end do
        end do
        end do


        CALL FACTOR_OF_TWO_TERMS(rhou_int,rho_int,u_cur)

        CALL FACTOR_OF_TWO_TERMS(rhov_int,rho_int,v_cur)

        CALL FACTOR_OF_TWO_TERMS(rhow_int,rho_int,w_cur)

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = ke * 0.5d0*rho_int(i,j,k)

        p_int(i,j,k) = (e_int(i,j,k) - ke) * (gama - 1.0D0)

        end do
        end do
        end do

        !k = NZF
        !do j = NYB,NYT
        !do i = NXL,NXR

        !rho_int(i,j,k)  = rho_int(i,j,NZB)

        !rhou_int(i,j,k) = rhou_int(i,j,NZB)

        !rhov_int(i,j,k) = rhov_int(i,j,NZB)

        !rhow_int(i,j,k) = rhow_int(i,j,NZB)

        !p_int(i,j,k) = p_int(i,j,NZB)

        !e_int(i,j,k) = e_int(i,j,NZB)

        !end do
        !end do

        !write(*,*) 'After Updating STEP1'
        !write(*,*) p_int(NXL,NYB,NZF),e_int(NXL,NYB,NZB),de1(NXL,NYB,NZB)


        !!LEFT BOTTOM
        !i = NXL

        !do j = NYB+1,NY1-1
        !do k = NZB,NZF-1

        !CALL Y_COORD(y,j)
        !y = abs(y)
        !v_in       = v_jet * (0.545d0 + 0.455d0*tanh((d-y)/(d/10.0d0)) )

        !rho_fac      = (1.0d0 + v_in/v_jet)*M*M*(gama-1.0d0)*0.5d0
        !rho_fac      = 0.0d0 - rho_fac*v_in/v_jet + 1.0d0
        !rho_fac      = 1.0d0/rho_fac
        !rho_in       = rho_fac * rho_atm

        !rho_int(i,j,k) = rho_int(i,j,k) *(1.0d0 - sig) + sig * rho_in

        !u_cur(i,j,k) = (1.0d0 - sig)*u_cur(i,j,k) + sig * v_in

        !rhou_int(i,j,k) = rho_int(i,j,k) * u_cur(i,j,k)

        !rhov_int(i,j,k) = rho_int(i,j,k) * v_cur(i,j,k)

        !rhow_int(i,j,k) = rho_int(i,j,k) * w_cur(i,j,k)

        !p_int(i,j,k) = p_int(i,j,k)*(1.0d0 - sig) + sig * p_atm

        !ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        !ke = 0.5d0 * ke * rho_int(i,j,k)

        !e_int(i,j,k) = p_int(i,j,k)/(gama - 1.0d0) + ke

        !end do
        !end do

        !i = NXL
        !k = NZF
        !do j = NYB+1,NY1-1

        !rho_int(i,j,k)  = rho_int(i,j,NZB)

        !rhou_int(i,j,k) = rhou_int(i,j,NZB)

        !rhov_int(i,j,k) = rhov_int(i,j,NZB)

        !rhow_int(i,j,k) = rhow_int(i,j,NZB)

        !p_int(i,j,k) = p_int(i,j,NZB)

        !e_int(i,j,k) = e_int(i,j,NZB)

        !end do


        !!LEFT TOP
        !i = NXL

        !do j = NY2+1,NYT-1
        !do k = NZB,NZF-1

        !CALL Y_COORD(y,j)
        !y = abs(y)
        !v_in       = v_jet * (0.545d0 + 0.455d0*tanh((d-y)/(d/10.0d0)) )

        !rho_fac      = (1.0d0 + v_in/v_jet)*M*M*(gama-1.0d0)*0.5d0
        !rho_fac      = 0.0d0 - rho_fac*v_in/v_jet + 1.0d0
        !rho_fac      = 1.0d0/rho_fac
        !rho_in       = rho_fac * rho_atm

        !rho_int(i,j,k) = rho_int(i,j,k) *(1.0d0 - sig) + sig * rho_in

        !u_cur(i,j,k) = (1.0d0 - sig)*u_cur(i,j,k) + sig * v_in

        !rhou_int(i,j,k) = rho_int(i,j,k) * u_cur(i,j,k)

        !rhov_int(i,j,k) = rho_int(i,j,k) * v_cur(i,j,k)

        !rhow_int(i,j,k) = rho_int(i,j,k) * w_cur(i,j,k)

        !p_int(i,j,k) = p_int(i,j,k)*(1.0d0 - sig) + sig * p_atm

        !ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        !ke = 0.5d0 * ke * rho_int(i,j,k)

        !e_int(i,j,k) = p_int(i,j,k)/(gama - 1.0d0) + ke

        !end do
        !end do

        !i = NXL
        !k = NZF
        !do j = NY2+1,NYT-1

        !rho_int(i,j,k)  = rho_int(i,j,NZB)

        !rhou_int(i,j,k) = rhou_int(i,j,NZB)

        !rhov_int(i,j,k) = rhov_int(i,j,NZB)

        !rhow_int(i,j,k) = rhow_int(i,j,NZB)

        !p_int(i,j,k) = p_int(i,j,NZB)

        !e_int(i,j,k) = e_int(i,j,NZB)

        !end do


        !INFLOW
        i = NXL

        do k = NZB,NZF
        do j = NYB,NYT

        !y          = (j - (NY1+NY2)/2) * dy1
        CALL Z_COORD(z,k)
        CALL Y_COORD(y,j)
        rad = sqrt(z*z+y*y)

        if (rad <= jet_r) then

        !v_in       = v_jet * (0.5d0 + 0.5d0*tanh((jet_r-rad)/(jet_r/10.0d0)) )
        !v_in       = v_jet * (0.5d0 + 0.5d0*tanh((jet_r-rad)/(jet_r/b_shear(j,k))) )
        v_in       = u_tar(j,k) !+ u_if(j,k)

        rho_fac      = 1.0d0/rho_in *(1.0d0 + (gama -1.0d0)*M*M*0.5d0) - 1.0d0/rho_c 
        rho_fac      = rho_fac*v_in/v_jet + 1.0d0/rho_c
        rho_fac      = rho_fac - (gama-1.0d0)*M*M*0.5d0*(v_in/v_jet)*(v_in/v_jet)/rho_in
        rho_inlet    = 1.0d0/rho_fac

        !rho_fac      = (1.0d0 + v_in/v_jet)*M*M*(gama-1.0d0)*0.5d0
        !rho_fac      = 0.0d0 - rho_fac*v_in/v_jet + 1.0d0
        !rho_fac      = 1.0d0/rho_fac
        !rho_in       = rho_fac * rho_atm
        !rho_in       = rho_atm !+ rho_if(j,k)

        rho_int(i,j,k) = rho_int(i,j,k) *(1.0d0 - sigi) + sigi * rho_inlet

        u_cur(i,j,k) = (1.0d0 - sigi)*u_cur(i,j,k) + sigi * v_in + u_if(j,k)

        rhou_int(i,j,k) = rho_int(i,j,k) * u_cur(i,j,k)

        v_cur(i,j,k) = (1.0d0 - sigi)*v_cur(i,j,k) + v_if(j,k) !+ sigi * v_if(j,k)

        rhov_int(i,j,k) = rho_int(i,j,k) * v_cur(i,j,k)

        w_cur(i,j,k) = (1.0d0 - sigi)*w_cur(i,j,k) + w_if(j,k)  !+ sigi * w_if(j,k)

        rhow_int(i,j,k) = rho_int(i,j,k) * w_cur(i,j,k)

        p_int(i,j,k) = p_int(i,j,k)*(1.0d0 - sigi) + sigi * p_in

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_int(i,j,k)

        e_int(i,j,k) = p_int(i,j,k)/(gama - 1.0d0) + ke

        else

        v_in       = v_jet * (0.5d0 + 0.5d0*tanh((jet_r-rad)/(jet_r/10.0d0)) )
        v_in       = v_in !+ u_if(j,k)

        rho_fac      = 1.0d0/rho_in *(1.0d0 + (gama -1.0d0)*M*M*0.5d0) - 1.0d0/rho_c 
        rho_fac      = rho_fac*v_in/v_jet + 1.0d0/rho_c
        rho_fac      = rho_fac - (gama-1.0d0)*M*M*0.5d0*(v_in/v_jet)*(v_in/v_jet)/rho_in
        rho_inlet    = 1.0d0/rho_fac

        rho_int(i,j,k) = rho_int(i,j,k) *(1.0d0 - sig) + sig * rho_inlet

        rhou_int(i,j,k) = rho_int(i,j,k) * u_cur(i,j,k)

        rhov_int(i,j,k) = rho_int(i,j,k) * v_cur(i,j,k)

        rhow_int(i,j,k) = rho_int(i,j,k) * w_cur(i,j,k)

        p_int(i,j,k) = p_int(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_int(i,j,k)

        e_int(i,j,k) = p_int(i,j,k)/(gama - 1.0d0) + ke

        end if

        end do
        end do

 


        !BOTTOM
        j = NYB

        do k = NZB,NZF
        do i = NXL+1,NXR-1

        rho_int(i,j,k) = rho_int(i,j,k) *(1.0d0 - sig) + sig * rho_c

        rhou_int(i,j,k) = rho_int(i,j,k) * u_cur(i,j,k)

        rhov_int(i,j,k) = rho_int(i,j,k) * v_cur(i,j,k)

        rhow_int(i,j,k) = rho_int(i,j,k) * w_cur(i,j,k)

        p_int(i,j,k) = p_int(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_int(i,j,k)

        e_int(i,j,k) = p_int(i,j,k)/(gama - 1.0d0) + ke

        end do
        end do

 


        !TOP
        j = NYT

        do k = NZB,NZF
        do i = NXL+1,NXR-1

        rho_int(i,j,k) = rho_int(i,j,k) *(1.0d0 - sig) + sig * rho_c

        rhou_int(i,j,k) = rho_int(i,j,k) * u_cur(i,j,k)

        rhov_int(i,j,k) = rho_int(i,j,k) * v_cur(i,j,k)

        rhow_int(i,j,k) = rho_int(i,j,k) * w_cur(i,j,k)

        p_int(i,j,k) = p_int(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_int(i,j,k)

        e_int(i,j,k) = p_int(i,j,k)/(gama - 1.0d0) + ke

        end do
        end do



        !BACK
        k = NZB

        do j = NYB+1,NYT-1
        do i = NXL+1,NXR-1

        rho_int(i,j,k) = rho_int(i,j,k) *(1.0d0 - sig) + sig * rho_c

        rhou_int(i,j,k) = rho_int(i,j,k) * u_cur(i,j,k)

        rhov_int(i,j,k) = rho_int(i,j,k) * v_cur(i,j,k)

        rhow_int(i,j,k) = rho_int(i,j,k) * w_cur(i,j,k)

        p_int(i,j,k) = p_int(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_int(i,j,k)

        e_int(i,j,k) = p_int(i,j,k)/(gama - 1.0d0) + ke

        end do
        end do

 


        !FRONT
        k = NZF

        do j = NYB+1,NYT-1
        do i = NXL+1,NXR-1

        rho_int(i,j,k) = rho_int(i,j,k) *(1.0d0 - sig) + sig * rho_c

        rhou_int(i,j,k) = rho_int(i,j,k) * u_cur(i,j,k)

        rhov_int(i,j,k) = rho_int(i,j,k) * v_cur(i,j,k)

        rhow_int(i,j,k) = rho_int(i,j,k) * w_cur(i,j,k)

        p_int(i,j,k) = p_int(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_int(i,j,k)

        e_int(i,j,k) = p_int(i,j,k)/(gama - 1.0d0) + ke

        end do
        end do


        CALL TEMPERATURE1

        !write(*,*) 'After Correction STEP 1'
        !write(*,*) p_int(NXL,NYB,NZF),e_int(NXL,NYB,NZF)

        END SUBROUTINE
        !_____________________________________________________________________

