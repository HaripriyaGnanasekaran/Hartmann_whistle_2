        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 16-02-2009        Last Modified on 16-05-2010

        !______________________________________________________________________





        !______________________________________________________________________
        MODULE y_momentum

        USE GRID_DIMENSIONS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: v_cur
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: v_int
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: v_new
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: v_mean
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: v_ms
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: v_tot
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: v_ts

        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: dv1
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: dv2
        !REAL (KIND=8), DIMENSION(1:NT,NXL:NXR), SAVE :: dev

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: v_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: v_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: v_z

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhov_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhov_int
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhov_new
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: rhov_int
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: rhov_new
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhov_tot

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: drhov1
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: drhov2

        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhov_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhov_y
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhov_z

        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: vrhov_xx
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: vrhov_yy
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: vrhov_zz

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vrhov_xx
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vrhov_yy
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vrhov_zz

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhovu_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhovv_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhovw_z

        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: v_xx
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: v_yy
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: v_zz

        REAL (KIND=8), DIMENSION(NYB:NYT,NZB:NZF)        , SAVE :: v_if


        END MODULE
        !______________________________________________________________________





