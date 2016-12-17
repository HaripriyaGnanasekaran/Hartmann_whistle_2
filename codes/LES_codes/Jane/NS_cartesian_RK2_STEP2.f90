        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 11-07-2009        Last Modified on 30-05-2010

        !_____________________________________________________________________




        !_____________________________________________________________________
        SUBROUTINE RK2_STEP2
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NY1,NY2
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dy1,dz1,dt,jet_r

        USE FLOW_PARAMETERS

        USE continuity, ONLY : rho_cur,rho_new,drho1,drho2,rho_if
        USE x_momentum
        USE y_momentum, ONLY : rhov_cur,rhov_new,drhov1,drhov2,v_cur
        USE y_momentum, ONLY : v_if
        USE z_momentum, ONLY : rhow_cur,rhow_new,drhow1,drhow2,w_cur
        USE z_momentum, ONLY : w_if
        USE pressure,   ONLY : p_new
        USE energy,     ONLY : e_cur,e_new,de1,de2
        IMPLICIT NONE

        INTEGER :: i,j,k

        REAL (KIND=8) :: ke,d,v_in,M,sig,sigi,u_in,rho_inlet
        REAL (KIND=8) :: x,y,z,rad,rho_fac

        c = sqrt(gama*p_in/rho_in)
        M = v_jet/c
        d = (NY2-NY1)*dy1*0.5d0

        sig     = 0.005d0
        sigi    = 0.05d0


        !ALLOCATE(rho_new(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(rhou_new(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(rhov_new(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(rhow_new(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(p_new(NXL:NXR,NYB:NYT,NZB:NZF))
        !ALLOCATE(e_new(NXL:NXR,NYB:NYT,NZB:NZF))

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF

        !rho_new(i,j,k)  = 0.0d0
        !rhou_new(i,j,k) = 0.0d0
        !rhov_new(i,j,k) = 0.0d0
        !rhow_new(i,j,k) = 0.0d0
        !p_new(i,j,k)    = 0.0d0
        !e_new(i,j,k)    = 0.0d0

        !end do
        !end do
        !end do

        !write(*,*) 'After Allocation STEP2' 
        !write(*,*) p_new(NXL,NYB,NZF),e_new(NXL,NYB,NZF)
        !!write(*,*) p_int(NXL,NYB,NZF),e_int(NXL,NYB,NZF)


        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        rho_new(i,j,k)  = rho_cur(i,j,k) - 0.5d0*(drho1(i,j,k) + drho2(i,j,k))*dt

        rhou_new(i,j,k) = rhou_cur(i,j,k) - 0.5d0*(drhou1(i,j,k) + drhou2(i,j,k))*dt

        rhov_new(i,j,k) = rhov_cur(i,j,k) - 0.5d0*(drhov1(i,j,k) + drhov2(i,j,k))*dt

        rhow_new(i,j,k) = rhow_cur(i,j,k) - 0.5d0*(drhow1(i,j,k) + drhow2(i,j,k))*dt

        e_new(i,j,k)    = e_cur(i,j,k) - 0.5d0*(de1(i,j,k)+de2(i,j,k))*dt 

        end do
        end do
        end do


        CALL FACTOR_OF_TWO_TERMS(rhou_new,rho_new,u_cur)

        CALL FACTOR_OF_TWO_TERMS(rhov_new,rho_new,v_cur)

        CALL FACTOR_OF_TWO_TERMS(rhow_new,rho_new,w_cur)

        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = ke*0.5d0 *rho_new(i,j,k)

        p_new(i,j,k) = (e_new(i,j,k) - ke) * (gama - 1.0D0)

        end do
        end do
        end do

        !k = NZF
        !do j = NYB,NYT
        !do i = NXL,NXR

        !rho_new(i,j,k)  = rho_new(i,j,NZB)

        !rhou_new(i,j,k) = rhou_new(i,j,NZB)

        !rhov_new(i,j,k) = rhov_new(i,j,NZB)

        !rhow_new(i,j,k) = rhow_new(i,j,NZB)

        !p_new(i,j,k)    = p_new(i,j,NZB)

        !e_new(i,j,k)    = e_new(i,j,NZB)

        !end do
        !end do

        !write(*,*) 'After Updating STEP2'
        !write(*,*) p_new(NXL,NYB,NZF),e_new(NXL,NYB,NZB),de1(NXL,NYB,NZB),de2(NXL,NYB,NZB)
        !write(*,*) rhou_new(NXL,NYB,NZF),drhou1(NXL,NYB,NZB),drhou2(NXL,NYB,NZB)

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

        !rho_new(i,j,k)  = rho_new(i,j,k) *(1.0d0 - sigi) + sigi * rho_in

        !u_cur(i,j,k) = (1.0d0 - sigi)*u_cur(i,j,k) + sig * v_in 

        !rhou_new(i,j,k) = rho_new(i,j,k) * u_cur(i,j,k)

        !rhov_new(i,j,k) = rho_new(i,j,k) * v_cur(i,j,k)

        !rhow_new(i,j,k) = rho_new(i,j,k) * w_cur(i,j,k)

        !p_new(i,j,k)    = p_new(i,j,k)*(1.0d0 - sig) + sig * p_atm

        !ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        !ke = 0.5d0 * ke * rho_new(i,j,k)

        !e_new(i,j,k) = p_new(i,j,k)/(gama - 1.0d0) + ke

        !end do
        !end do

        !i = NXL
        !k = NZF
        !do j = NYB+1,NY1-1

        !rho_new(i,j,k)  = rho_new(i,j,NZB)

        !rhou_new(i,j,k) = rhou_new(i,j,NZB)

        !rhov_new(i,j,k) = rhov_new(i,j,NZB)

        !rhow_new(i,j,k) = rhow_new(i,j,NZB)

        !p_new(i,j,k)    = p_new(i,j,NZB)

        !e_new(i,j,k)    = e_new(i,j,NZB)

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

        !rho_new(i,j,k)  = rho_new(i,j,k) *(1.0d0 - sig) + sig * rho_in

        !u_cur(i,j,k) = (1.0d0 - sig)*u_cur(i,j,k) + sig * v_in 

        !rhou_new(i,j,k) = rho_new(i,j,k) * u_cur(i,j,k)

        !rhov_new(i,j,k) = rho_new(i,j,k) * v_cur(i,j,k)

        !rhow_new(i,j,k) = rho_new(i,j,k) * w_cur(i,j,k)

        !p_new(i,j,k)    = p_new(i,j,k)*(1.0d0 - sig) + sig * p_atm

        !ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        !ke = 0.5d0 * ke * rho_new(i,j,k)

        !e_new(i,j,k) = p_new(i,j,k)/(gama - 1.0d0) + ke

        !end do
        !end do

        !i = NXL
        !k = NZF
        !do j = NY2+1,NYT-1

        !rho_new(i,j,k)  = rho_new(i,j,NZB)

        !rhou_new(i,j,k) = rhou_new(i,j,NZB)

        !rhov_new(i,j,k) = rhov_new(i,j,NZB)

        !rhow_new(i,j,k) = rhow_new(i,j,NZB)

        !p_new(i,j,k)    = p_new(i,j,NZB)

        !e_new(i,j,k)    = e_new(i,j,NZB)

        !end do


        !INFLOW
        i = NXL

        do k = NZB,NZF
        do j = NYB,NYT

        !y          = (j - (NY1+ NY2)/2) * dy1
        CALL Z_COORD(z,k)
        CALL Y_COORD(y,j)
        rad = sqrt(z*z+y*y)

        if ( rad <= jet_r ) then

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
        !rho_in       = rho_atm  !+ rho_if(j,k)

        rho_new(i,j,k) = rho_new(i,j,k) *(1.0d0 - sigi) + sigi * rho_inlet

        u_cur(i,j,k) = (1.0d0 - sigi)*u_cur(i,j,k) + sigi * v_in + u_if(j,k) 

        rhou_new(i,j,k) = rho_new(i,j,k) * u_cur(i,j,k)

        v_cur(i,j,k) = (1.0d0 - sigi)*v_cur(i,j,k) + v_if(j,k) !+ sigi * v_if(j,k)

        rhov_new(i,j,k) = rho_new(i,j,k) * v_cur(i,j,k)

        w_cur(i,j,k) = (1.0d0 - sigi)*w_cur(i,j,k) + w_if(j,k) !+ sigi * w_if(j,k)

        rhow_new(i,j,k) = rho_new(i,j,k) * w_cur(i,j,k)

        p_new(i,j,k) = p_new(i,j,k)*(1.0d0 - sigi) + sigi * p_in

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_new(i,j,k)

        e_new(i,j,k) = p_new(i,j,k)/(gama - 1.0d0) + ke

        else

        v_in       = v_jet * (0.5d0 + 0.5d0*tanh((jet_r-rad)/(jet_r/10.0d0)) )
        v_in       = v_in !+ u_if(j,k)

        rho_fac      = 1.0d0/rho_in *(1.0d0 + (gama -1.0d0)*M*M*0.5d0) - 1.0d0/rho_c 
        rho_fac      = rho_fac*v_in/v_jet + 1.0d0/rho_c
        rho_fac      = rho_fac - (gama-1.0d0)*M*M*0.5d0*(v_in/v_jet)*(v_in/v_jet)/rho_in
        rho_inlet    = 1.0d0/rho_fac

        rho_new(i,j,k) = rho_new(i,j,k) *(1.0d0 - sig) + sig * rho_inlet

        rhou_new(i,j,k) = rho_new(i,j,k) * u_cur(i,j,k)

        rhov_new(i,j,k) = rho_new(i,j,k) * v_cur(i,j,k)

        rhow_new(i,j,k) = rho_new(i,j,k) * w_cur(i,j,k)

        p_new(i,j,k) = p_new(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_new(i,j,k)

        e_new(i,j,k) = p_new(i,j,k)/(gama - 1.0d0) + ke

        end if

        end do
        end do


        
        
        !BOTTOM
        j = NYB

        do k = NZB,NZF
        do i = NXL+1,NXR-1

        rho_new(i,j,k)  = rho_new(i,j,k) *(1.0d0 - sig) + sig * rho_c

        rhou_new(i,j,k) = rho_new(i,j,k) * u_cur(i,j,k)

        rhov_new(i,j,k) = rho_new(i,j,k) * v_cur(i,j,k)

        rhow_new(i,j,k) = rho_new(i,j,k) * w_cur(i,j,k)

        p_new(i,j,k)    = p_new(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_new(i,j,k)

        e_new(i,j,k) = p_new(i,j,k)/(gama - 1.0d0) + ke

        end do
        end do

 


        !TOP
        j = NYT

        do k = NZB,NZF
        do i = NXL+1,NXR-1

        rho_new(i,j,k)  = rho_new(i,j,k) *(1.0d0 - sig) + sig * rho_c

        rhou_new(i,j,k) = rho_new(i,j,k) * u_cur(i,j,k)

        rhov_new(i,j,k) = rho_new(i,j,k) * v_cur(i,j,k)

        rhow_new(i,j,k) = rho_new(i,j,k) * w_cur(i,j,k)

        p_new(i,j,k)    = p_new(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_new(i,j,k)

        e_new(i,j,k) = p_new(i,j,k)/(gama - 1.0d0) + ke

        end do
        end do



        !BACK
        k = NZB

        do j = NYB+1,NYT-1
        do i = NXL+1,NXR-1

        rho_new(i,j,k)  = rho_new(i,j,k) *(1.0d0 - sig) + sig * rho_c

        rhou_new(i,j,k) = rho_new(i,j,k) * u_cur(i,j,k)

        rhov_new(i,j,k) = rho_new(i,j,k) * v_cur(i,j,k)

        rhow_new(i,j,k) = rho_new(i,j,k) * w_cur(i,j,k)

        p_new(i,j,k)    = p_new(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_new(i,j,k)

        e_new(i,j,k) = p_new(i,j,k)/(gama - 1.0d0) + ke

        end do
        end do

 


        !FRONT
        k = NZF

        do j = NYB+1,NYT-1
        do i = NXL+1,NXR-1

        rho_new(i,j,k)  = rho_new(i,j,k) *(1.0d0 - sig) + sig * rho_c

        rhou_new(i,j,k) = rho_new(i,j,k) * u_cur(i,j,k)

        rhov_new(i,j,k) = rho_new(i,j,k) * v_cur(i,j,k)

        rhow_new(i,j,k) = rho_new(i,j,k) * w_cur(i,j,k)

        p_new(i,j,k)    = p_new(i,j,k)*(1.0d0 - sig) + sig * p_c

        ke = u_cur(i,j,k)**2 + v_cur(i,j,k)**2 + w_cur(i,j,k)**2
        ke = 0.5d0 * ke * rho_new(i,j,k)

        e_new(i,j,k) = p_new(i,j,k)/(gama - 1.0d0) + ke

        end do
        end do


        !write(*,*) u_cur(NXL,NY1,NY1),u_if(NY1,NY1)

        !write(*,*) 'After Correction STEP 2'
        !write(*,*) p_new(NXL,NYB,NZB),e_new(NXL,NYB,NZB)
        !write(*,*) u_cur(NXL,NYB,NZB),rho_new(NXL,NYB,NZB),rhou_new(NXL,NYB,NZB)
        !write(*,*) p_new(NXL,NYB,NZF),e_new(NXL,NYB,NZF)
        !write(*,*) u_cur(NXL,NYB,NZF),rho_new(NXL,NYB,NZF),rhou_new(NXL,NYB,NZF)


        END SUBROUTINE
        !_____________________________________________________________________

