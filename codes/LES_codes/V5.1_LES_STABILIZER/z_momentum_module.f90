        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 16-02-2009        Last Modified on 10-05-2010

        !______________________________________________________________________





        !______________________________________________________________________
        MODULE z_momentum

        USE GRID_DIMENSIONS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: w_cur
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: w_int
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: w_new
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: w_mean
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: w_ms
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: w_tot
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: w_ts

        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: dw1
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: dw2
        !REAL (KIND=8), DIMENSION(1:NT,NXL:NXR), SAVE :: dev

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: w_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: w_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: w_z

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhow_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhow_int
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhow_new
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: rhow_int
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: rhow_new
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhow_tot

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: drhow1
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: drhow2

        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhow_x
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhow_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhow_z

        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: vrhow_xx
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: vrhow_yy
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: vrhow_zz

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vrhow_xx
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vrhow_yy
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vrhow_zz

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhowu_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhowv_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhoww_z

        REAL (KIND=8), DIMENSION(NYB:NYT,NZB:NZF)        , SAVE :: w_if

        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: w_xx
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: w_yy
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: w_zz

        END MODULE
        !______________________________________________________________________





