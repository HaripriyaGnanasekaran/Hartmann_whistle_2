        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 26-06-2012        Last Modified on 26-06-2012

        !______________________________________________________________________







        !______________________________________________________________________
        SUBROUTINE CALCULATE_LIGHTHILL1

        USE GRID_DIMENSIONS
        USE GRID_SIZE_TIME_STEP
        USE FLOW_PARAMETERS

        USE continuity
        USE x_momentum
        USE y_momentum
        USE z_momentum
        USE pressure

        USE lighthill_source

        IMPLICIT NONE

        CHARACTER(len = 40) :: fileoutput
        REAL (KIND=8) :: x,y,z,c_o,c_o_sqr
        REAL (KIND=8) :: omega

        INTEGER :: i,j,k

        c_o_sqr = gama*p_c/rho_c
        c_o     = sqrt(c_o_sqr)
  

!!!     X Derivatives      
!!$OMP PARALLEL SHARED(rhouu_x,rhouv_y,rhouw_z,p_x,rho_x,lh_s_x) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s_x(i,j,k) =  rhouu_x(i,j,k)  + rhouv_y(i,j,k)  + rhouw_z(i,j,k)
        lh_s_x(i,j,k) =  lh_s_x(i,j,k) - vrhou_xx(i,j,k) - vrhou_yy(i,j,k) - vrhou_zz(i,j,k)
        !lh_s_x(i,j,k) =  lh_s_x(i,j,k) + p_x(i,j,k) - c_o_sqr * (rho_x(i,j,k) - rho_mean_x(i,j,k))
        lh_s_x(i,j,k) =  lh_s_x(i,j,k) + p_x(i,j,k)- c_o_sqr * rho_x(i,j,k) 
        end do
        end do
        end do     
!!$OMP END DO
!!$OMP END PARALLEL

        !i = (NXL+NXR)/2
        !j = NY1
        !k = NY1
        !write(*,*) rhouu_x(i,j,k) + rhouv_y(i,j,k) + rhouw_z(i,j,k)
        !write(*,*) p_x(i,j,k) 
        !write(*,*) c_o_sqr*rho_x(i,j,k),c_o_sqr*rho_mean_x(i,j,k)
        !write(*,*) c_o_sqr*(rho_x(i,j,k) - rho_mean_x(i,j,k))


!!!     Y Derivatives      
!!$OMP PARALLEL SHARED(rhovu_x,rhovv_y,rhovw_z,p_y,rho_y,lh_s_y) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s_y(i,j,k) =  rhovu_x(i,j,k) + rhovv_y(i,j,k) + rhovw_z(i,j,k)
        lh_s_y(i,j,k) =  lh_s_y(i,j,k) - vrhov_xx(i,j,k) - vrhov_yy(i,j,k) - vrhov_zz(i,j,k)
        !lh_s_y(i,j,k) =  lh_s_y(i,j,k) + p_y(i,j,k) - c_o_sqr * (rho_y(i,j,k) - rho_mean_y(i,j,k))
        lh_s_y(i,j,k) =  lh_s_y(i,j,k) + p_y(i,j,k) - c_o_sqr * rho_y(i,j,k) 
        end do
        end do
        end do     
!!$OMP END DO
!!$OMP END PARALLEL





!!!     Z Derivatives    
!!$OMP PARALLEL SHARED(rhowu_x,rhowv_y,rhoww_z,p_z,rho_z,lh_s_z) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s_z(i,j,k) =  rhowu_x(i,j,k) + rhowv_y(i,j,k) + rhoww_z(i,j,k)
        lh_s_z(i,j,k) =  lh_s_z(i,j,k) - vrhow_xx(i,j,k) - vrhow_yy(i,j,k) - vrhow_zz(i,j,k)
        !lh_s_z(i,j,k) =  lh_s_z(i,j,k) + p_z(i,j,k) - c_o_sqr * (rho_z(i,j,k) - rho_mean_z(i,j,k))
        lh_s_z(i,j,k) =  lh_s_z(i,j,k) + p_z(i,j,k) - c_o_sqr * rho_z(i,j,k) 
        end do
        end do
        end do     
!!$OMP END DO
!!$OMP END PARALLEL





        END  SUBROUTINE
        !______________________________________________________________________








        !______________________________________________________________________
        SUBROUTINE CALCULATE_LIGHTHILL2

        USE GRID_DIMENSIONS
        USE GRID_SIZE_TIME_STEP
        USE FLOW_PARAMETERS

        USE continuity
        USE x_momentum
        USE y_momentum
        USE z_momentum
        USE pressure

        USE lighthill_source

        IMPLICIT NONE

        CHARACTER(len = 40) :: fileoutput
        REAL (KIND=8) :: x,y,z,c_o,c_o_sqr
        REAL (KIND=8) :: omega

        INTEGER :: i,j,k

        c_o_sqr = gama*p_c/rho_c
        c_o     = sqrt(c_o_sqr)

        !write(*,*) c_o_sqr,c_o


!!!     X Derivatives      
!!$OMP PARALLEL SHARED(rhouu_x,rhouv_y,rhouw_z,p_x,rho_x,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s_x(i,j,k) =  lh_s_x(i,j,k) + rhouu_x(i,j,k)  + rhouv_y(i,j,k)  + rhouw_z(i,j,k)
        lh_s_x(i,j,k) =  lh_s_x(i,j,k) - vrhou_xx(i,j,k) - vrhou_yy(i,j,k) - vrhou_zz(i,j,k)
        !lh_s_x(i,j,k) =  lh_s_x(i,j,k) + p_x(i,j,k) - c_o_sqr * (rho_x(i,j,k) - rho_mean_x(i,j,k))
        lh_s_x(i,j,k) =  lh_s_x(i,j,k) + p_x(i,j,k) - c_o_sqr * rho_x(i,j,k) 
        end do
        end do
        end do     
!!$OMP END DO
!!$OMP END PARALLEL


        !i = (NXL+NXR)/2
        !j = NY1
        !k = NY1
        !write(*,*) rhouu_x(i,j,k) + rhouv_y(i,j,k) + rhouw_z(i,j,k)
        !write(*,*) p_x(i,j,k) 
        !write(*,*) c_o_sqr*rho_x(i,j,k), c_o_sqr*rho_mean_x(i,j,k)



!!!     Y Derivatives      
!!$OMP PARALLEL SHARED(rhovu_x,rhovv_y,rhovw_z,p_y,rho_y,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s_y(i,j,k) =  lh_s_y(i,j,k) + rhovu_x(i,j,k) + rhovv_y(i,j,k) + rhovw_z(i,j,k)
        lh_s_y(i,j,k) =  lh_s_y(i,j,k) - vrhov_xx(i,j,k) - vrhov_yy(i,j,k) - vrhov_zz(i,j,k)
        !lh_s_y(i,j,k) =  lh_s_y(i,j,k) + p_y(i,j,k) - c_o_sqr*(rho_y(i,j,k) - rho_mean_y(i,j,k))
        lh_s_y(i,j,k) =  lh_s_y(i,j,k) + p_y(i,j,k)  - c_o_sqr * rho_y(i,j,k) 
        end do
        end do
        end do     
!!$OMP END DO
!!$OMP END PARALLEL






!!!     Z Derivatives    
!!$OMP PARALLEL SHARED(rhowu_x,rhowv_y,rhoww_z,p_z,rho_z,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s_z(i,j,k) =  lh_s_z(i,j,k) + rhowu_x(i,j,k) + rhowv_y(i,j,k) + rhoww_z(i,j,k)
        lh_s_z(i,j,k) =  lh_s_z(i,j,k) - vrhow_xx(i,j,k) - vrhow_yy(i,j,k) - vrhow_zz(i,j,k)
        !lh_s_z(i,j,k) =  lh_s_z(i,j,k) + p_z(i,j,k) - c_o_sqr*(rho_z(i,j,k) - rho_mean_z(i,j,k))
        lh_s_z(i,j,k) =  lh_s_z(i,j,k) + p_z(i,j,k) - c_o_sqr * rho_z(i,j,k) 
        end do
        end do
        end do     
!!$OMP END DO
!!$OMP END PARALLEL



        END  SUBROUTINE
        !______________________________________________________________________








        !______________________________________________________________________
        SUBROUTINE CALCULATE_LIGHTHILL(l)

        USE GRID_DIMENSIONS
        USE GRID_SIZE_TIME_STEP
        USE FLOW_PARAMETERS

        USE continuity
        USE x_momentum
        USE y_momentum
        USE z_momentum
        USE pressure

        USE lighthill_source

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 40) :: fileoutput
        REAL (KIND=8) :: x,y,z,c_o
        REAL (KIND=8) :: omega

        INTEGER :: i,j,k

        c_o = sqrt(gama*p_c/rho_c)


!!!     X Derivatives      
!!$OMP PARALLEL SHARED(lh_s_f,lh_s_b,lh_s_x) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s_f(i,j,k) =  lh_s_x(i,j,k)
        lh_s_b(i,j,k) =  lh_s_x(i,j,k)
        end do
        end do
        end do     
!!$OMP END DO
!!$OMP END PARALLEL


        CALL HIXON_BACKDIF_X(lh_s_b)
        CALL HIXON_FORWARD_X(lh_s_f)


!!$OMP PARALLEL SHARED(lh_s,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s(i,j,k) = lh_s_b(i,j,k)+ lh_s_f(i,j,k)
        end do
        end do
        end do 
!!$OMP END DO
!!$OMP END PARALLEL





!!!     Y Derivatives      
!!$OMP PARALLEL SHARED(lh_s_f,lh_s_b,lh_s_y) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s_f(i,j,k) =  lh_s_y(i,j,k) 
        lh_s_b(i,j,k) =  lh_s_y(i,j,k)
        end do
        end do
        end do     
!!$OMP END DO
!!$OMP END PARALLEL


        CALL HIXON_BACKDIF_Y(lh_s_b)
        CALL HIXON_FORWARD_Y(lh_s_f)

!!$OMP PARALLEL SHARED(lh_s,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s(i,j,k) = lh_s(i,j,k) + lh_s_b(i,j,k)+ lh_s_f(i,j,k)
        end do
        end do
        end do 
!!$OMP END DO
!!$OMP END PARALLEL






!!!     Z Derivatives    
!!$OMP PARALLEL SHARED(lh_s_f,lh_s_b,lh_s_z) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s_f(i,j,k) =  lh_s_z(i,j,k) 
        lh_s_b(i,j,k) =  lh_s_z(i,j,k)
        end do
        end do
        end do     
!!$OMP END DO
!!$OMP END PARALLEL


        CALL HIXON_BACKDIF_Z(lh_s_b)
        CALL HIXON_FORWARD_Z(lh_s_f)

!!$OMP PARALLEL SHARED(lh_s,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        lh_s(i,j,k) = lh_s(i,j,k) + lh_s_b(i,j,k)+ lh_s_f(i,j,k)
        lh_s(i,j,k) = 0.25d0 * lh_s(i,j,k)
        end do
        end do
        end do 
!!$OMP END DO
!!$OMP END PARALLEL


        END  SUBROUTINE
        !______________________________________________________________________








        !______________________________________________________________________
        SUBROUTINE CALCULATE_MEAN_P_RHO_GRAD(l)

        USE GRID_DIMENSIONS
        USE GRID_SIZE_TIME_STEP
        USE FLOW_PARAMETERS

        USE continuity
        USE pressure

        USE lighthill_source

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 40) :: fileoutput
        REAL (KIND=8) :: x,y,z,c_o
        REAL (KIND=8) :: omega

        INTEGER :: i,j,k

        c_o = sqrt(gama*p_c/rho_c)


        !!!!! MEAN PRESSURE GRAD


        CALL COPY(p_tot,lh_s_b)
        CALL COPY(lh_s_b,lh_s_f)

        CALL HIXON_BACKDIF_X(lh_s_b)
        CALL HIXON_FORWARD_X(lh_s_f)


!!$OMP PARALLEL SHARED(lh_s,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        p_mean_x(i,j,k) = 0.5d0*(lh_s_b(i,j,k)+ lh_s_f(i,j,k))
        end do
        end do
        end do 
!!$OMP END DO
!!$OMP END PARALLEL


        CALL COPY(p_tot,lh_s_b)
        CALL COPY(lh_s_b,lh_s_f)

        CALL HIXON_BACKDIF_Y(lh_s_b)
        CALL HIXON_FORWARD_Y(lh_s_f)


!!$OMP PARALLEL SHARED(lh_s,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        p_mean_y(i,j,k) = 0.5d0*(lh_s_b(i,j,k)+ lh_s_f(i,j,k))
        end do
        end do
        end do 
!!$OMP END DO
!!$OMP END PARALLEL

        CALL COPY(p_tot,lh_s_b)
        CALL COPY(lh_s_b,lh_s_f)

        CALL HIXON_BACKDIF_Z(lh_s_b)
        CALL HIXON_FORWARD_Z(lh_s_f)


!!$OMP PARALLEL SHARED(lh_s,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        p_mean_z(i,j,k) = 0.5d0*(lh_s_b(i,j,k)+ lh_s_f(i,j,k))
        end do
        end do
        end do 
!!$OMP END DO
!!$OMP END PARALLEL







        !!!!! MEAN DENSITY GRAD


        CALL COPY(rho_tot,lh_s_b)
        CALL COPY(lh_s_b,lh_s_f)

        CALL HIXON_BACKDIF_X(lh_s_b)
        CALL HIXON_FORWARD_X(lh_s_f)


!!$OMP PARALLEL SHARED(lh_s,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        rho_mean_x(i,j,k) = 0.5d0*(lh_s_b(i,j,k)+ lh_s_f(i,j,k))
        end do
        end do
        end do 
!!$OMP END DO
!!$OMP END PARALLEL


        CALL COPY(rho_tot,lh_s_b)
        CALL COPY(lh_s_b,lh_s_f)

        CALL HIXON_BACKDIF_Y(lh_s_b)
        CALL HIXON_FORWARD_Y(lh_s_f)


!!$OMP PARALLEL SHARED(lh_s,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        rho_mean_y(i,j,k) = 0.5d0*(lh_s_b(i,j,k)+ lh_s_f(i,j,k))
        end do
        end do
        end do 
!!$OMP END DO
!!$OMP END PARALLEL

        CALL COPY(rho_tot,lh_s_b)
        CALL COPY(lh_s_b,lh_s_f)

        CALL HIXON_BACKDIF_Z(lh_s_b)
        CALL HIXON_FORWARD_Z(lh_s_f)


!!$OMP PARALLEL SHARED(lh_s,lh_s_f,lh_s_b) PRIVATE(i,j,k)
!!$OMP DO SCHEDULE(DYNAMIC)
        do k = NZB,NZF
        do j = NYB,NYT
        do i = NXL,NXR
        rho_mean_z(i,j,k) = 0.5d0*(lh_s_b(i,j,k)+ lh_s_f(i,j,k))
        end do
        end do
        end do 
!!$OMP END DO
!!$OMP END PARALLEL



        END  SUBROUTINE
        !______________________________________________________________________






