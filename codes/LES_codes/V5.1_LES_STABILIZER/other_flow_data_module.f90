        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 16-02-2009        Last Modified on 10-05-2010

        !______________________________________________________________________




        MODULE METRICS
        USE GRID_DIMENSIONS
        IMPLICIT NONE
        REAL(KIND=8), DIMENSION(NXL:NXR), SAVE :: dx_deta
        REAL(KIND=8), DIMENSION(NYB:NYT), SAVE :: dy_deta
        REAL(KIND=8), DIMENSION(NZB:NZF), SAVE :: dz_deta
        END MODULE
        !____________________________________________________________________




        !______________________________________________________________________
        MODULE pressure
        USE GRID_DIMENSIONS
        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_int
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_new
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: p_int
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: p_new
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_mean
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_ms
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_tot
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_ts

        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: dp1
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: dp2

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_z

        END MODULE
        !______________________________________________________________________





        
        !______________________________________________________________________
        MODULE vorticity
        USE GRID_DIMENSIONS
        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),SAVE :: omegax
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),SAVE :: omegay
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),SAVE :: omegaz

        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT) :: omega_xyprof

        END MODULE
        !______________________________________________________________________




        
        !______________________________________________________________________
        MODULE sound_field
        USE GRID_DIMENSIONS
        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF),SAVE :: spl_f

        END MODULE
        !______________________________________________________________________





        
        !______________________________________________________________________
        MODULE kirchoff
        USE GRID_DIMENSIONS
        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NYB:NYT,NZB:NZF),SAVE :: p_px      !Pressure Positive normal y surface
        REAL (KIND=8), DIMENSION(NYB:NYT,NZB:NZF),SAVE :: p_nx      !Pressure Negative normal y surface
        REAL (KIND=8), DIMENSION(NYB:NYT,NZB:NZF),SAVE :: dpdn_px     !Pressure Normal gradience Positive normal y surface
        REAL (KIND=8), DIMENSION(NYB:NYT,NZB:NZF),SAVE :: dpdn_nx     !Pressure Normal gradience Negative normal y surface

        REAL (KIND=8), DIMENSION(NXL:NXR,NZB:NZF),SAVE :: p_py      !Pressure Positive normal y surface
        REAL (KIND=8), DIMENSION(NXL:NXR,NZB:NZF),SAVE :: p_ny      !Pressure Negative normal y surface
        REAL (KIND=8), DIMENSION(NXL:NXR,NZB:NZF),SAVE :: dpdn_py     !Pressure Normal gradience Positive normal y surface
        REAL (KIND=8), DIMENSION(NXL:NXR,NZB:NZF),SAVE :: dpdn_ny     !Pressure Normal gradience Negative normal y surface

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT),SAVE :: p_pz      !Pressure Positive normal z surface
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT),SAVE :: p_nz      !Pressure Negative normal z surface
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT),SAVE :: dpdn_pz     !Pressure Normal gradience Positive normal z surface
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT),SAVE :: dpdn_nz     !Pressure Normal gradience Negative normal z surface


        END MODULE
        !______________________________________________________________________



        !______________________________________________________________________
        MODULE cross_momentum

        USE GRID_DIMENSIONS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: uv_tot
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: uw_tot
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vw_tot

         REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: turb_ke_tot

        END MODULE
        !______________________________________________________________________

 
      


        !______________________________________________________________________
        MODULE lighthill_source

        USE GRID_DIMENSIONS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: lh_s
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: lh_s_tot
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: lh_s_ts


        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: lh_s_f
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: lh_s_b

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: lh_s_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: lh_s_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: lh_s_z

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_mean_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_mean_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: p_mean_z

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rho_mean_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rho_mean_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rho_mean_z

        END MODULE
        !______________________________________________________________________



 
        !______________________________________________________________________
        MODULE Monitor_Points
        USE GRID_DIMENSIONS
        IMPLICIT NONE
        INTEGER, PARAMETER :: NP_MP = 11
        INTEGER, PARAMETER :: NT_MP = 200

        INTEGER, DIMENSION(1:NP_MP) :: NX_MP = (/ 64,96,215,96,96,96,96,215,215,215,215 /)
        INTEGER, DIMENSION(1:NP_MP) :: NY_MP = (/ 0 ,0  ,0 ,-8,-8,8 ,8 ,- 8,- 8, 8 ,8  /)
        INTEGER, DIMENSION(1:NP_MP) :: NZ_MP = (/ 0 ,0  ,0 ,-8,8 ,-8,8 ,- 8, 8 ,-8, 8  /) 

        REAL (KIND=8), DIMENSION(1:NP_MP,1:NT_MP), SAVE :: rho_mp
        REAL (KIND=8), DIMENSION(1:NP_MP,1:NT_MP), SAVE :: u_mp
        REAL (KIND=8), DIMENSION(1:NP_MP,1:NT_MP), SAVE :: v_mp
        REAL (KIND=8), DIMENSION(1:NP_MP,1:NT_MP), SAVE :: w_mp
        REAL (KIND=8), DIMENSION(1:NP_MP,1:NT_MP), SAVE :: p_mp

        END MODULE
        !______________________________________________________________________

        

       
        
        !______________________________________________________________________
        MODULE Monitor_Lines
        USE GRID_DIMENSIONS
        IMPLICIT NONE
        INTEGER, PARAMETER :: NL_ML = 5
        INTEGER, PARAMETER :: NT_ML = 2000

        !INTEGER, DIMENSION(1:NP_MP) :: NX_MP = (/ 64,96,215,96,96,96,96,215,215,215,215 /)
        !INTEGER, DIMENSION(1:NL_ML) :: NY_ML = (/ 0,-8,-8, 8,8 /)
        !INTEGER, DIMENSION(1:NL_ML) :: NZ_ML = (/ 0,-8, 8,-8,8 /) 
        INTEGER, DIMENSION(1:NL_ML) :: NY_ML = (/ 0, 0,0,-8,8 /)
        INTEGER, DIMENSION(1:NL_ML) :: NZ_ML = (/ 0,-8,8, 0,0 /) 

        REAL (KIND=8), DIMENSION(1:NL_ML,NXL:NXR,1:NT_ML), SAVE :: rho_ml
        REAL (KIND=8), DIMENSION(1:NL_ML,NXL:NXR,1:NT_ML), SAVE :: u_ml
        REAL (KIND=8), DIMENSION(1:NL_ML,NXL:NXR,1:NT_ML), SAVE :: v_ml
        REAL (KIND=8), DIMENSION(1:NL_ML,NXL:NXR,1:NT_ML), SAVE :: w_ml
        REAL (KIND=8), DIMENSION(1:NL_ML,NXL:NXR,1:NT_ML), SAVE :: p_ml
        REAL (KIND=8), DIMENSION(1:NL_ML,NXL:NXR,1:NT_ML), SAVE :: lh_s_ml

        END MODULE
        !______________________________________________________________________

        


        
!        !______________________________________________________________________
!        MODULE turbulence_module
!        USE GRID_DIMENSIONS
!        IMPLICIT NONE!

!        REAL (KIND=8),                SAVE :: A    	!X direction
!        REAL (KIND=8),                SAVE :: A_R  	!R direction
!        REAL (KIND=8),                SAVE :: A_THETA	!THETA direction
!        REAL (KIND=8),                SAVE :: A_shear	!THETA direction
!        REAL (KIND=8), DIMENSION(1:8),SAVE :: alpha
!        REAL (KIND=8), DIMENSION(1:8),SAVE :: psi

!        END MODULE
!        !______________________________________________________________________


        
        !______________________________________________________________________
        MODULE turbulence_freund
        USE GRID_DIMENSIONS
        IMPLICIT NONE
        INTEGER, PARAMETER :: NTF = 4  !Half number of modes
	INTEGER, PARAMETER :: NTA = 2  !Number of azimuthal modes NTF*NTA=total number of modes

        REAL (KIND=8), DIMENSION(1:NTF,0:NTA-1),SAVE :: A    	
        REAL (KIND=8),                SAVE :: dA  	
	REAL (KIND=8),                SAVE :: d_alpha
        REAL (KIND=8), DIMENSION(1:NTF,0:NTA-1),SAVE :: alpha
        REAL (KIND=8), DIMENSION(1:NTF,0:NTA-1),SAVE :: phi
        REAL (KIND=8), DIMENSION(1:NTF,0:NTA-1),SAVE :: psi

        END MODULE
        !______________________________________________________________________


