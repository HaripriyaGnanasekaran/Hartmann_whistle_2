        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 13-07-2009        Last Modified on 16-02-2010

        !_____________________________________________________________________


        !______________________________________________________________________
        SUBROUTINE RK2(ts)

        USE GRID_DIMENSIONS,   ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE continuity,    ONLY : rho_cur,rho_new,rho_int
        USE x_momentum,    ONLY : u_cur
        USE y_momentum,    ONLY : v_cur
        USE z_momentum,    ONLY : w_cur
        USE pressure,      ONLY : p_cur,p_new,p_int
        USE energy,        ONLY : e_cur,e_new,e_int
        USE x_momentum,    ONLY : rhou_cur,rhou_new,rhou_int
        USE x_momentum,    ONLY : vrhou_xx,vrhou_yy,vrhou_zz
        USE y_momentum,    ONLY : rhov_cur,rhov_new,rhov_int
        USE y_momentum,    ONLY : vrhov_xx,vrhov_yy,vrhov_zz
        USE z_momentum,    ONLY : rhow_cur,rhow_new,rhow_int
        USE z_momentum,    ONLY : vrhow_xx,vrhow_yy,vrhow_zz
       
        IMPLICIT NONE       

        INTEGER, INTENT(INOUT) :: ts

        REAL (KIND=8) :: x,y,z
        
        INTEGER :: i,j,k
        INTEGER, DIMENSION(2) :: ch
      
        !write(*,*) 'RK2 Step 1'  
       
        if (MOD(ts,4) == 0) then
        ch(1) = 0
        ch(2) = 1
        end if 
        if (MOD(ts,4) == 1) then
        ch(1) = 1
        ch(2) = 0
        end if
        if (MOD(ts,4) == 2) then
        ch(1) = 0
        ch(2) = 1
        end if
        if (MOD(ts,4) == 3) then
        ch(1) = 1
        ch(2) = 0
        end if

        CALL EVALUATION1(ch)

        CALL BOUNDARY_EDGE1(ch)

        CALL BOUNDARY_CORNER1(ch)

        CALL BOUNDARY1(ts)

        CALL RK2_STEP1

        CALL VORTICITY1

        CALL CALCULATE_KIRCHOFF1

        CALL CALCULATE_LIGHTHILL1


        !DEALLOCATE(vrhou_xx)
        !DEALLOCATE(vrhov_xx)
        !DEALLOCATE(vrhow_xx)
        !DEALLOCATE(vrhou_yy)
        !DEALLOCATE(vrhov_yy)
        !DEALLOCATE(vrhow_yy)
        !DEALLOCATE(vrhou_zz)
        !DEALLOCATE(vrhov_zz)
        !DEALLOCATE(vrhow_zz)


        CALL FACTOR_OF_TWO_TERMS(rhou_int,rho_int,u_cur)

        CALL FACTOR_OF_TWO_TERMS(rhov_int,rho_int,v_cur)

        CALL FACTOR_OF_TWO_TERMS(rhow_int,rho_int,w_cur)



        !write(*,*) 'STEP 1 Complete '

        if (MOD(ts,4)== 0) then
        ch(1) = 1
        ch(2) = 0
        else if (MOD(ts,4) == 1) then
        ch(1) = 0
        ch(2) = 1
        else if (MOD(ts,4) == 2) then
        ch(1) = 1
        ch(2) = 0
        else
        ch(1) = 0
        ch(2) = 1
        end if

        !write(*,*) ' Begin Differentiation'
        CALL EVALUATION2(ch)
        !write(*,*) 'Step 2 Derivatives calculated'

        CALL BOUNDARY_EDGE2(ch)

        CALL BOUNDARY_CORNER2(ch)

        CALL BOUNDARY2(ts)

        CALL RK2_STEP2

        CALL VORTICITY2

        CALL CALCULATE_LIGHTHILL2


        !DEALLOCATE(vrhou_xx)
        !DEALLOCATE(vrhov_xx)
        !DEALLOCATE(vrhow_xx)
        !DEALLOCATE(vrhou_yy)
        !DEALLOCATE(vrhov_yy)
        !DEALLOCATE(vrhow_yy)
        !DEALLOCATE(vrhou_zz)
        !DEALLOCATE(vrhov_zz)
        !DEALLOCATE(vrhow_zz)

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

        !CALL COPY(rho_int,rho_new)
        !CALL COPY(rhou_int,rhou_new)
        !CALL COPY(rhov_int,rhov_new)
        !CALL COPY(rhow_int,rhow_new)
        !CALL COPY(p_int,p_new)
        !CALL COPY(e_int,e_new)


        !DEALLOCATE(rho_int)
        !DEALLOCATE(rhou_int)
        !DEALLOCATE(rhov_int)
        !DEALLOCATE(rhow_int)
        !DEALLOCATE(p_int)
        !DEALLOCATE(e_int)

        if (ch(1) == 0) then
        CALL FILTERING_X(rho_new)
        CALL FILTERING_X(rhou_new)
        CALL FILTERING_X(rhov_new)
        CALL FILTERING_X(rhow_new)
        CALL FILTERING_X(p_new)
        CALL FILTERING_X(e_new)

        CALL FILTERING_Y(rho_new)
        CALL FILTERING_Y(rhou_new)
        CALL FILTERING_Y(rhov_new)
        CALL FILTERING_Y(rhow_new)
        CALL FILTERING_Y(p_new)
        CALL FILTERING_Y(e_new)

        CALL FILTERING_Z(rho_new)
        CALL FILTERING_Z(rhou_new)
        CALL FILTERING_Z(rhov_new)
        CALL FILTERING_Z(rhow_new)
        CALL FILTERING_Z(p_new)
        CALL FILTERING_Z(e_new)

        else
        CALL FILTERING_X(rho_new)
        CALL FILTERING_X(rhou_new)
        CALL FILTERING_X(rhov_new)
        CALL FILTERING_X(rhow_new)
        CALL FILTERING_X(p_new)
        CALL FILTERING_X(e_new)

        CALL FILTERING_Z(rho_new)
        CALL FILTERING_Z(rhou_new)
        CALL FILTERING_Z(rhov_new)
        CALL FILTERING_Z(rhow_new)
        CALL FILTERING_Z(p_new)
        CALL FILTERING_Z(e_new)


        CALL FILTERING_Y(rho_new)
        CALL FILTERING_Y(rhou_new)
        CALL FILTERING_Y(rhov_new)
        CALL FILTERING_Y(rhow_new)
        CALL FILTERING_Y(p_new)
        CALL FILTERING_Y(e_new)

 

 
        end if
        
        
        CALL BUFFER_FILTER(rho_new)
        CALL BUFFER_FILTER(rhou_new)
        CALL BUFFER_FILTER(rhov_new)
        CALL BUFFER_FILTER(rhow_new)
        CALL BUFFER_FILTER(e_new)
        CALL BUFFER_FILTER(p_new)

        CALL FACTOR_OF_TWO_TERMS(rhou_new,rho_new,u_cur)

        CALL FACTOR_OF_TWO_TERMS(rhov_new,rho_new,v_cur)

        CALL FACTOR_OF_TWO_TERMS(rhow_new,rho_new,w_cur)

        CALL TEMPERATURE2

        !write(*,*) 'STEP 2 Complete '
       
        CALL COPY(rho_new,rho_cur)
        CALL COPY(rhou_new,rhou_cur)
        CALL COPY(rhov_new,rhov_cur)
        CALL COPY(rhow_new,rhow_cur)
        CALL COPY(p_new,p_cur)
        CALL COPY(e_new,e_cur)

        CALL CALCULATE_KIRCHOFF2

        !DEALLOCATE(rho_new)
        !DEALLOCATE(rhou_new)
        !DEALLOCATE(rhov_new)
        !DEALLOCATE(rhow_new)
        !DEALLOCATE(p_new)
        !DEALLOCATE(e_new)

        END SUBROUTINE
        !______________________________________________________________________



