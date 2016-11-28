        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 16-02-2010        Last Modified on 10-05-2010

        !______________________________________________________________________





        !______________________________________________________________________
        MODULE energy
        USE GRID_DIMENSIONS
        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_int
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_new

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), sAVE :: T

        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: e_int
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: e_new

        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_mean
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_ms
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_tot
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_ts

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: de1
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: de2

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: T_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: T_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: T_z  

        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_x
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_y
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: e_z  
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: efl_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: efl_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: efl_z

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vke_xx
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vke_yy
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vke_zz

        END MODULE
        !______________________________________________________________________




        
 
