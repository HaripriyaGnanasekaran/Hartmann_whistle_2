        ! Purpose:
        ! This is the main code to solve the problem of 3D
        ! Round jets using full NS equations 

        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 15-08-2008        Last Modified on 03-07-2015

        !_____________________________________________________________________


        !Main Program
        !______________________________________________________________________
        PROGRAM Acoustic_Scattering

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_SIZE_TIME_STEP, ONLY : dt,jet_d
        USE FLOW_PARAMETERS
        USE x_momentum, ONLY : u_cur,rhou_cur !,du1,du2
        USE y_momentum, ONLY : v_cur,rhov_cur !,dv1,dv2
        USE z_momentum, ONLY : w_cur,rhow_cur !,dw1,dw2
        USE continuity, ONLY : rho_cur !,drho1,drho2
        USE pressure,   ONLY : p_cur !,dp1,dp2
        USE energy,     ONLY : e_cur
        !USE vorticity,  ONLY : omega


        IMPLICIT NONE
        !INTEGER, PARAMETER :: NP = 8000
        INTEGER (KIND=4) :: i,j,k,mi,meanorig,meanstart,laststep,temp
        INTEGER (KIND=4) :: development_flowthrough,no_of_flowthrough_mean,tft
        INTEGER (KIND=4) :: data_write_every,mean_write_every
        INTEGER :: l,m
        INTEGER :: rem 
        REAL (KIND=8) :: t,xr
        REAL(KIND=8), DIMENSION(NXL:NXR-1):: xr_arr
        REAL :: start_time,end_time
        WRITE(*,*) ''
        WRITE(*,*) ''
        WRITE(*,*)'******************************************************************'
        WRITE(*,*)'                  LARGE EDDY SIMULATION                           '
        WRITE(*,*)'******************************************************************'
        WRITE(*,*)''
        WRITE(*,*)''
        WRITE(*,*)'SIMULATING ROUND JET (COMPRESSIBLE)'
        WRITE(*,*)'MACH NUMBER = ',v_jet/sqrt(gama*R*T_in)
        WRITE(*,*)'REYNOLS NO. = ',rho_in*v_jet*jet_d/mu
        i = 0
        j = 0
        k = 1
        mi = 0

        OPEN(UNIT=18,FILE='input.dat',ACTION='READ')
        READ(18,*) i
        READ(18,*) j
        READ(18,*) development_flowthrough
        READ(18,*) no_of_flowthrough_mean
        READ(18,*) data_write_every,mean_write_every
        CLOSE(18)

        CALL DEVELOPMENTMODULE(development_flowthrough,tft)
        meanstart=tft
        CALL FLOWTHROUGHMODULE(no_of_flowthrough_mean,tft)
        laststep=meanstart+tft
        meanorig = meanstart       
        WRITE(*,*) 'SIMULATION BEGINS AT         =',I,' steps'
        IF(j.ne.0)THEN
        WRITE(*,*) 'MEAN DATA IS READ FROM       =',j,' steps'
        END IF
        WRITE(*,*) 'DATA COLLECTION STARTS AT    =',meanstart,' steps'
        write(*,*) 'SIMULATION ENDS AT TIME STEP =',laststep,' steps'
        WRITE(*,*)'******************************************************************'
        

        write(*,*)  'Initialization'
        CALL INITIALIZATION
        CALL DENSITY_INITIALIZATION
        CALL U_VELOCITY_INITIALIZATION
        CALL V_VELOCITY_INITIALIZATION
        CALL PRESSURE_INITIALIZATION
        CALL ENERGY_INITIALIZATION
        CALL VORTICITY_INITIALIZATION
        CALL LIGHTHILL_INITIALIZATION
        CALL INITIALIZE_FILTERING_COEFF_X
        CALL INITIALIZE_FILTERING_COEFF_Y
        CALL INITIALIZE_FILTERING_COEFF_Z
        CALL GRID_INITIALIZATION
        CALL INITIALIZE_BUFFER_FILTER
        CALL METRICS_X
        CALL METRICS_Y
        CALL METRICS_Z
        CALL INITIALIZE_RANDOM_NUMBERS_FREUND
        !CALL DE_INITIALIZATION

        write(*,*) 'Initialization Done'

        temp=i

	IF (i .ne. 0) then
        CALL READ_FILE(i)
        CALL READ_DENSITY(i)
        CALL READ_U_VELOCITY(i)
        CALL READ_V_VELOCITY(i)
        CALL READ_W_VELOCITY(i)
        CALL READ_PRESSURE(i)
        CALL READ_ENERGY(i)
        CALL READ_VORTICITY(i)
        CALL READ_LIGHTHILL_SOURCE(i)
	CALL READ_INFLOW_FREUND(i)
        END IF
        
        CALL PRODUCT_OF_TWO_TERMS(rho_cur,u_cur,rhou_cur)
        CALL PRODUCT_OF_TWO_TERMS(rho_cur,v_cur,rhov_cur)
        CALL PRODUCT_OF_TWO_TERMS(rho_cur,w_cur,rhow_cur)
        
        IF(j.ne.0)THEN
	i = j
        CALL READ_MEAN_DENSITY(i)
        CALL READ_MEAN_U_VELOCITY(i)
        CALL READ_MEAN_V_VELOCITY(i)
        CALL READ_MEAN_W_VELOCITY(i)
        CALL READ_MEAN_PRESSURE(i)
        CALL READ_MEAN_ENERGY(i)
        CALL READ_MEAN_LIGHTHILL_SOURCE(i)
        CALL READ_MEAN_CROSS_MOMENTUM(i)
        CALL READ_MEAN_RHOU_MOMENTUM(i)
        CALL READ_MEAN_RHOV_MOMENTUM(i)
        CALL READ_MEAN_RHOW_MOMENTUM(i)
        CALL READ_MEAN_TURB_KE(i)
        !CALL CALCULATE_MEAN_P_RHO_GRAD(l)
        CALL READ_TOTAL(i)
	END IF

	i =temp
        CALL CPU_TIME(start_time)
        do
        if(i >= laststep) exit 
        !SOLUTION BLOWUP CONTROLLER
        !**************************
        if(i .gt. meanorig)then
        meanstart = meanorig
        else
        meanstart=0
        end if
        !***************************
        CALL TOTAL_VALUE
        CALL TOTAL_SQUARE_VALUE
        CALL INFLOW_FREUND(i)
        CALL RK2(i)

        i = i + 1
        mi = mi + 1
       

        IF(i.le.meanstart)THEN
        rem = MOD(i,data_write_every)
        ELSE
        rem = MOD(i,mean_write_every)
        END IF

        if (rem == 0 ) then

        CALL WRITE_INFLOW_FREUND(i)

        CALL WRITE_DENSITY(i)
        CALL WRITE_U_VELOCITY(i)
        CALL WRITE_V_VELOCITY(i)
        CALL WRITE_W_VELOCITY(i)
        CALL WRITE_PRESSURE(i)
        CALL WRITE_ENERGY(i)
        CALL WRITE_VORTICITY(i)
        CALL WRITE_LIGHTHILL_SOURCE(i)
        
        CALL WRITE_DENSITY_XY_PROFILE(i)
        CALL WRITE_U_VELOCITY_XY_PROFILE(i)
        CALL WRITE_V_VELOCITY_XY_PROFILE(i)
        CALL WRITE_W_VELOCITY_XY_PROFILE(i)
        CALL WRITE_PRESSURE_XY_PROFILE(i)
        CALL WRITE_ENERGY_XY_PROFILE(i)
        CALL WRITE_LIGHTHILL_SOURCE_XY_PROFILE(i)
        CALL WRITE_TEMPERATURE_XY_PROFILE(i)
        CALL WRITE_VORTICITY_XY_PROFILE(i)

        CALL WRITE_DENSITY_XZ_PROFILE(i)
        CALL WRITE_U_VELOCITY_XZ_PROFILE(i)
        CALL WRITE_V_VELOCITY_XZ_PROFILE(i)
        CALL WRITE_W_VELOCITY_XZ_PROFILE(i)
        CALL WRITE_PRESSURE_XZ_PROFILE(i)
        CALL WRITE_ENERGY_XZ_PROFILE(i)
        CALL WRITE_VORTICITY_XZ_PROFILE(i)
        
        CALL WRITE_DENSITY_YZ_PROFILE(i)
        CALL WRITE_U_VELOCITY_YZ_PROFILE(i)
        CALL WRITE_V_VELOCITY_YZ_PROFILE(i)
        CALL WRITE_W_VELOCITY_YZ_PROFILE(i)
        CALL WRITE_PRESSURE_YZ_PROFILE(i)
        CALL WRITE_ENERGY_YZ_PROFILE(i)
        CALL WRITE_VORTICITY_YZ_PROFILE(i)
	CALL WRITE_INFLOW_FREUND(i)

        CALL CALCULATE_MONITORPOINTS(mi)
        CALL CALCULATE_MONITORLINES(mi)
        CALL CALCULATE_LIGHTHILL(i)

        WRITE(*,*) 'STEP',i,'     COMPLETED...   FILES ARE WRITTEN'
        end if
        
        OPEN(UNIT=19,FILE='run_monitor_log.dat',ACTION='WRITE',POSiTION='APPEND')
        write(19,*) 'Step',i,'u(0,0,0)',u_cur(0,0,0)
        CLOSE(19)
         

	if(i > meanstart) then
        j = i - meanstart
       
                if(i .lt. meanorig)then
                j = 4000
                write(*,*) 'success'
                end if

        rem = MOD(j,mean_write_every)
       
        if (rem == 0 ) then
        !
        CALL MEAN_VALUE(j)
        CALL MEAN_SQUARE_VALUE(j)
        !
        CALL WRITE_MEAN_DENSITY(j)
        CALL WRITE_MEAN_U_VELOCITY(j)
        CALL WRITE_MEAN_V_VELOCITY(j)
        CALL WRITE_MEAN_W_VELOCITY(j)
        CALL WRITE_MEAN_PRESSURE(j)
        CALL WRITE_MEAN_ENERGY(j)
        CALL WRITE_MEAN_LIGHTHILL_SOURCE(j)
        CALL WRITE_MEAN_CROSS_MOMENTUM(j)
        CALL WRITE_MEAN_RHOU_MOMENTUM(j)
        CALL WRITE_MEAN_RHOV_MOMENTUM(j)
        CALL WRITE_MEAN_RHOW_MOMENTUM(j)
        CALL WRITE_MEAN_TURB_KE(j)
        !
        CALL READ_TOTAL(j)
        !
        CALL WRITE_MONITORPOINTS(j)
        CALL WRITE_MONITORLINES(j)
        WRITE(*,*) 'MEAN DATA WRITTEN FOR',J,'SAMPLES'
        !
        mi = 0
        end if
	end if


        end do


        CALL CPU_TIME(end_time)
        PRINT *, 'TIME', end_time - start_time
        WRITE(*,*) 'LES SIMULATION - ROUNDJET: COMPLETED SUCCESSFULLY'
	END PROGRAM Acoustic_Scattering






        !______________________________________________________________________
        SUBROUTINE VALUE (m, c)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: m
        CHARACTER, INTENT(OUT) :: c

        SELECT CASE (m)
        CASE (0)
                c = '0'
        CASE (1) 
                c = '1'
        CASE (2)
                c = '2'
        CASE (3)
                c = '3'
        CASE (4)
                c = '4'
        CASE (5)
                c = '5'
        CASE (6)
                c = '6'
        CASE (7)
                c = '7'
        CASE (8)
                c = '8'
        CASE (9) 
                c = '9'

        END SELECT
        c = c
        END SUBROUTINE
        !______________________________________________________________________


        !______________________________________________________________________
        SUBROUTINE COPY(val_cur,dval)
        
        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),    INTENT(IN)  :: val_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),    INTENT(OUT) :: dval
       
        INTEGER :: i,j,k
        
!$OMP PARALLEL SHARED(val_cur,dval) PRIVATE(i,j,k)

!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        dval(i,j,k) =  val_cur(i,j,k)
        end do
        end do
        end do     
!$OMP END DO
!$OMP END PARALLEL

        END SUBROUTINE
        !______________________________________________________________________

