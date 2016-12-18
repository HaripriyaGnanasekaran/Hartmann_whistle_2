        PROGRAM Jane_Watsons ! Just ANother Efficient and Wide Algorithm To Solve Navier Stokes

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_SIZE_TIME_STEP, ONLY : dt,jet_d
        USE FLOW_PARAMETERS
        USE x_momentum, ONLY : u_cur,rhou_cur !,du1,du2
        USE y_momentum, ONLY : v_cur,rhov_cur !,dv1,dv2
        USE z_momentum, ONLY : w_cur,rhow_cur !,dw1,dw2
        USE continuity, ONLY : rho_cur !,drho1,drho2
        USE pressure,   ONLY : p_cur !,dp1,dp2
        USE energy,     ONLY : e_cur
        


        IMPLICIT NONE
        
        INTEGER (KIND=4) :: i,j,k,mi,meanorig,meanstart,laststep,temp
        INTEGER (KIND=4) :: development_flowthrough,no_of_flowthrough_mean,tft
        INTEGER (KIND=4) :: data_write_every,mean_write_every
        INTEGER :: l,m
        INTEGER :: rem 
        REAL (KIND=8) :: t,xr
        REAL(KIND=8), DIMENSION(NXL:NXR-1):: xr_arr
        REAL :: start_time,end_time
        
        WRITE(*,*)''
        WRITE(*,*)''
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
        WRITE(*,*) 'SIMULATION ENDS AT TIME STEP =',laststep,' steps'
        WRITE(*,*)'******************************************************************'
        

       
		CALL INIT
        temp=i
		IF (i .ne. 0) then
			CALL READFILE(i)
        END IF
      
        CALL PRODUCT_OF_TWO_TERMS(rho_cur,u_cur,rhou_cur)
        CALL PRODUCT_OF_TWO_TERMS(rho_cur,v_cur,rhov_cur)
        CALL PRODUCT_OF_TWO_TERMS(rho_cur,w_cur,rhow_cur)

		IF(j.ne.0)THEN
			i = j
			CALL READMEAN(i)
		END IF
		i =temp
        CALL CPU_TIME(start_time)

        DO
        	IF(i >= 151) EXIT

        		!SOLUTION BLOWUP CONTROLLER
        		!**************************
        
				IF(i .GT. meanorig)THEN
        		meanstart = meanorig
        		ELSE
        		meanstart=0
        		END IF
        		!***************************
        
			CALL TOTAL_VALUE
        	CALL TOTAL_SQUARE_VALUE
        	CALL INFLOW_FREUND(i)
        	CALL RK2(i)
        	i = i + 1
        	mi = mi + 1
			
			if (i>50) then
			call gen_data_for_lcs(i-49)
			endif
			
	        IF(i.le.meanstart)THEN
        		rem = MOD(i,data_write_every)
        	ELSE
        		rem = MOD(i,mean_write_every)
        	END IF

        	IF (rem == 0 ) THEN
				!CALL WRITEFILE(i)
	        END IF
        
			OPEN(UNIT=19,FILE='run_monitor_log.dat',ACTION='WRITE',POSiTION='APPEND')
			write(19,*) 'Step',i,'u(0,0,0)',u_cur(0,0,0)
			CLOSE(19)
         

			IF(i > meanstart) THEN
        		j = i - meanstart
       
                IF(i .lt. meanorig)THEN
                	j = 4000
                	!write(*,*) 'success'
                END IF

        		rem = MOD(j,mean_write_every)
       
        		IF (rem == 0 ) THEN
					!CALL WRITEMEAN(j)
			        mi = 0
        		END IF
        	END IF

        END DO


        CALL CPU_TIME(end_time)
        PRINT *, 'TIME', end_time - start_time
        WRITE(*,*) 'LES SIMULATION - ROUNDJET: COMPLETED SUCCESSFULLY'
        END PROGRAM Jane_Watsons





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



		SUBROUTINE INIT
		
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
			
		END SUBROUTINE INIT
		
		
		SUBROUTINE READFILE(k)
			integer, intent(in) :: k
			CALL READ_FILE(k)
        	CALL READ_DENSITY(k)
        	CALL READ_U_VELOCITY(k)
        	CALL READ_V_VELOCITY(k)
        	CALL READ_W_VELOCITY(k)
        	CALL READ_PRESSURE(k)
        	CALL READ_ENERGY(k)
    	    CALL READ_VORTICITY(k)
	        CALL READ_LIGHTHILL_SOURCE(k)
			CALL READ_INFLOW_FREUND(k)
	
		END SUBROUTINE READFILE
	
	
		SUBROUTINE READMEAN(k)
			integer, intent(in) :: k
			CALL READ_MEAN_DENSITY(k)
      		CALL READ_MEAN_U_VELOCITY(k)
        	CALL READ_MEAN_V_VELOCITY(k)
        	CALL READ_MEAN_W_VELOCITY(k)
        	CALL READ_MEAN_PRESSURE(k)
        	CALL READ_MEAN_ENERGY(k)
        	CALL READ_MEAN_LIGHTHILL_SOURCE(k)
        	CALL READ_MEAN_CROSS_MOMENTUM(k)
        	CALL READ_MEAN_RHOU_MOMENTUM(k)
        	CALL READ_MEAN_RHOV_MOMENTUM(k)
        	CALL READ_MEAN_RHOW_MOMENTUM(k)
        	CALL READ_MEAN_TURB_KE(k)
        	CALL READ_TOTAL(k)
		
		END SUBROUTINE READMEAN
		
		SUBROUTINE WRITEFILE(k)
			integer, intent(in) :: k
			CALL WRITE_INFLOW_FREUND(k)
			CALL WRITE_DENSITY(k)
			CALL WRITE_U_VELOCITY(k)
			CALL WRITE_V_VELOCITY(k)
			CALL WRITE_W_VELOCITY(k)
			CALL WRITE_PRESSURE(k)
			CALL WRITE_ENERGY(k)
			CALL WRITE_VORTICITY(k)
			CALL WRITE_LIGHTHILL_SOURCE(k)

			CALL WRITE_DENSITY_XY_PROFILE(k)
			CALL WRITE_U_VELOCITY_XY_PROFILE(k)
			CALL WRITE_V_VELOCITY_XY_PROFILE(k)
			CALL WRITE_W_VELOCITY_XY_PROFILE(k)
			CALL WRITE_PRESSURE_XY_PROFILE(k)
			CALL WRITE_ENERGY_XY_PROFILE(k)
			CALL WRITE_LIGHTHILL_SOURCE_XY_PROFILE(k)
			CALL WRITE_TEMPERATURE_XY_PROFILE(k)
			CALL WRITE_VORTICITY_XY_PROFILE(k)

			CALL WRITE_DENSITY_XZ_PROFILE(k)
			CALL WRITE_U_VELOCITY_XZ_PROFILE(k)
			CALL WRITE_V_VELOCITY_XZ_PROFILE(k)
			CALL WRITE_W_VELOCITY_XZ_PROFILE(k)
			CALL WRITE_PRESSURE_XZ_PROFILE(k)
			CALL WRITE_ENERGY_XZ_PROFILE(k)
			CALL WRITE_VORTICITY_XZ_PROFILE(k)

			CALL WRITE_DENSITY_YZ_PROFILE(k)
			CALL WRITE_U_VELOCITY_YZ_PROFILE(k)
			CALL WRITE_V_VELOCITY_YZ_PROFILE(k)
			CALL WRITE_W_VELOCITY_YZ_PROFILE(k)
			CALL WRITE_PRESSURE_YZ_PROFILE(k)
			CALL WRITE_ENERGY_YZ_PROFILE(k)
			CALL WRITE_VORTICITY_YZ_PROFILE(k)
			CALL WRITE_INFLOW_FREUND(k)

			CALL CALCULATE_MONITORPOINTS(mi)
			CALL CALCULATE_MONITORLINES(mi)
			CALL CALCULATE_LIGHTHILL(k)

			WRITE(*,*) 'STEP',k,'     COMPLETED...   FILES ARE WRITTEN'
		
		END SUBROUTINE WRITEFILE
		
		
		SUBROUTINE WRITEMEAN(k)
			integer, intent(in) :: k
			CALL MEAN_VALUE(k)
			CALL MEAN_SQUARE_VALUE(k)
			!
			CALL WRITE_MEAN_DENSITY(k)
			CALL WRITE_MEAN_U_VELOCITY(k)
			CALL WRITE_MEAN_V_VELOCITY(k)
			CALL WRITE_MEAN_W_VELOCITY(k)
			CALL WRITE_MEAN_PRESSURE(k)
			CALL WRITE_MEAN_ENERGY(k)
			CALL WRITE_MEAN_LIGHTHILL_SOURCE(k)
			CALL WRITE_MEAN_CROSS_MOMENTUM(k)
			CALL WRITE_MEAN_RHOU_MOMENTUM(k)
			CALL WRITE_MEAN_RHOV_MOMENTUM(k)
			CALL WRITE_MEAN_RHOW_MOMENTUM(k)
			CALL WRITE_MEAN_TURB_KE(k)
			!
			CALL READ_TOTAL(k)
			!
			CALL WRITE_MONITORPOINTS(k)
			CALL WRITE_MONITORLINES(k)
			WRITE(*,*) 'MEAN DATA WRITTEN FOR',k,'SAMPLES'
		
		END SUBROUTINE WRITEMEAN
