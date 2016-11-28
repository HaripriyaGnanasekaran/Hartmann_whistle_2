        ! This program contains the modules which contains the
        ! coefficeints for filtering
        
        !______________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 13-07-2009        Last Modified on 26-05-2010

        !______________________________________________________________________


        !______________________________________________________________________
        MODULE FILTERING_PARAMETERS

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        IMPLICIT NONE
        SAVE

        REAL (KIND=8), DIMENSION(NXL:NXR-4,1:3) :: A_coef_x
        REAL (KIND=8), DIMENSION(NXL:NXR-5)     :: inv_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT-4,1:3) :: A_coef_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT-5)     :: inv_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NZB:NZF-4,1:3) :: A_coef_z
        REAL (KIND=8), DIMENSION(NXL:NXR,NZB:NZF-5)     :: inv_z 

        !REAL (KIND=8), DIMENSION(NZB::NZF-1,1:3) :: A1_coef_z
        !REAL (KIND=8), DIMENSION(NZB::NZF-2)     :: inv1_z
        !REAL (KIND=8), DIMENSION(NZB::NZF-1)     :: B_1
        !REAL (KIND=8), DIMENSION(NZB::NZF-1,1:3) :: A2_coef_z
        !REAL (KIND=8), DIMENSION(NZB::NZF-2)     :: inv2_z
        !REAL (KIND=8), DIMENSION(NZB::NZF-1,1:2) :: B_2

        REAL (KIND=8), PARAMETER :: alphax = 0.475d0
        REAL (KIND=8), DIMENSION(NXL:NXR) :: alphay != 0.498d0 !0.475d0 !0.498d0
        REAL (KIND=8), DIMENSION(NXL:NXR) :: alphaz != 0.498d0 !0.475d0 !0.498d0
        
        REAL (KIND=8), PARAMETER :: ax = (5.0d0 + 6.0d0*alphax)/8.0d0
        REAL (KIND=8), PARAMETER :: bx = (1.0d0 + 2.0d0*alphax)/2.0d0*0.5d0
        REAL (KIND=8), PARAMETER :: cx = 0.0d0 - (1.0d0 - 2.0d0*alphax)/8.0d0*0.5d0

        REAL (KIND=8), DIMENSION(NXL:NXR) :: ay != (5.0d0 + 6.0d0*alphay)/8.0d0
        REAL (KIND=8), DIMENSION(NXL:NXR) :: by != (1.0d0 + 2.0d0*alphay)/2.0d0*0.5d0
        REAL (KIND=8), DIMENSION(NXL:NXR) :: cy != 0.0d0 - (1.0d0 - 2.0d0*alphay)/8.0d0*0.5d0

        REAL (KIND=8), DIMENSION(NXL:NXR) :: az != (5.0d0 + 6.0d0*alphaz)/8.0d0
        REAL (KIND=8), DIMENSION(NXL:NXR) :: bz != (1.0d0 + 2.0d0*alphaz)/2.0d0*0.5d0
        REAL (KIND=8), DIMENSION(NXL:NXR) :: cz != 0.0d0 - (1.0d0 - 2.0d0*alphaz)/8.0d0*0.5d0


        END MODULE FILTERING_PARAMETERS
        !______________________________________________________________________


        !______________________________________________________________________
        MODULE BUFFER_FILTER_PARAMETERS

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NXB
        IMPLICIT NONE
        SAVE

        REAL (KIND=8), DIMENSION(NXB+1:NXR-1) :: x_buf


        END MODULE BUFFER_FILTER_PARAMETERS
        !______________________________________________________________________



