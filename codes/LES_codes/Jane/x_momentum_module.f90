        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 16-02-2010        Last Modified on 16-05-2010

        !______________________________________________________________________





        !______________________________________________________________________
        MODULE x_momentum

        USE GRID_DIMENSIONS

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: u_cur
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: u_int
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: u_new
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: u_mean    
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: u_ms
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: u_tot
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: u_ts

        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: du1
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: du2
        !REAL (KIND=8), DIMENSION(1:NT,NXL:NXR), SAVE :: deu

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: u_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: u_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: u_z

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhou_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhou_int
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhou_new
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: rhou_int
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: rhou_new
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhou_tot

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: drhou1
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: drhou2

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhou_x
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhou_y
        !REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhou_z

        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: vrhou_xx
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: vrhou_yy
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: vrhou_zz

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vrhou_xx
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vrhou_yy
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: vrhou_zz

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhouu_x
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhouv_y
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE :: rhouw_z

        REAL (KIND=8), DIMENSION(NYB:NYT,NZB:NZF)        , SAVE :: u_if
        REAL (KIND=8), DIMENSION(NYB:NYT,NZB:NZF)        , SAVE :: u_tar
        REAL (KIND=8), DIMENSION(NYB:NYT,NZB:NZF)        , SAVE :: b_shear

        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: u_xx
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: u_yy
        !REAL (KIND=8), ALLOCATABLE, DIMENSION(:,:,:), SAVE :: u_zz

        END MODULE
        !______________________________________________________________________


