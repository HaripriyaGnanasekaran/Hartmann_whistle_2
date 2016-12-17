        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 02-02-2010        Last Modified on 17-05-2010

        !_____________________________________________________________________






        !______________________________________________________________________
        SUBROUTINE WRITE_DENSITY_YZ_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity,          ONLY : rho_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 24) :: fileoutput
        REAL (KIND=8) :: x,y,z
        INTEGER :: m
        INTEGER :: UN
        INTEGER :: TE
        INTEGER :: HU
        INTEGER :: TH
        INTEGER :: TT
        INTEGER :: LC
        INTEGER :: i,j,k

        m    = 0
        UN   = 0
        TE   = 0
        HU   = 0
        TH   = 0
        TT   = 0
        LC   = 0

        fileoutput = 'den_yzprof_000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(17:17))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(16:16))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(15:15))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(14:14))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(13:13))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(12:12))

        OPEN(UNIT = 61, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        i = (NXL+NXR)/2

        do k = NZB,NZF
        do j = NYB,NYT

        write(61,*)rho_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 61)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_U_VELOCITY_YZ_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE x_momentum,          ONLY : u_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 24) :: fileoutput
        REAL (KIND=8) :: x,y,z
        INTEGER :: m
        INTEGER :: UN
        INTEGER :: TE
        INTEGER :: HU
        INTEGER :: TH
        INTEGER :: TT
        INTEGER :: LC
        INTEGER :: i,j,k

        m    = 0
        UN   = 0
        TE   = 0
        HU   = 0
        TH   = 0
        TT   = 0
        LC   = 0

        fileoutput = 'uvel_yzprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(17:17))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(16:16))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(15:15))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(14:14))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(13:13))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(12:12))

        OPEN(UNIT = 62, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        i = (NXL+NXR)/2
        
        do k = NZB,NZF
        do j = NYB,NYT

        write(62,*)u_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 62)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_V_VELOCITY_YZ_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE y_momentum,          ONLY : v_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 24) :: fileoutput
        REAL (KIND=8) :: x,y,z
        INTEGER :: m
        INTEGER :: UN
        INTEGER :: TE
        INTEGER :: HU
        INTEGER :: TH
        INTEGER :: TT
        INTEGER :: LC
        INTEGER :: i,j,k

        m    = 0
        UN   = 0
        TE   = 0
        HU   = 0
        TH   = 0
        TT   = 0
        LC   = 0

        fileoutput = 'vvel_yzprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(17:17))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(16:16))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(15:15))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(14:14))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(13:13))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(12:12))

        OPEN(UNIT = 63, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        i = (NXL+NXR)/2
        
        do k = NZB,NZF
        do j = NYB,NYT

        write(63,*)v_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 63)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_W_VELOCITY_YZ_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE z_momentum,          ONLY : w_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 24) :: fileoutput
        REAL (KIND=8) :: x,y,z
        INTEGER :: m
        INTEGER :: UN
        INTEGER :: TE
        INTEGER :: HU
        INTEGER :: TH
        INTEGER :: TT
        INTEGER :: LC
        INTEGER :: i,j,k

        m    = 0
        UN   = 0
        TE   = 0
        HU   = 0
        TH   = 0
        TT   = 0
        LC   = 0

        fileoutput = 'wvel_yzprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(17:17))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(16:16))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(15:15))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(14:14))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(13:13))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(12:12))

        OPEN(UNIT = 64, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        i = (NXL+NXR)/2
        
        do k = NZB,NZF
        do j = NYB,NYT

        write(64,*)w_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 64)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_PRESSURE_YZ_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE pressure,            ONLY : p_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 24) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
        REAL (KIND=8) :: x
        REAL (KIND=8) :: y    
        INTEGER :: m
        INTEGER :: UN
        INTEGER :: TE
        INTEGER :: HU
        INTEGER :: TH
        INTEGER :: TT
        INTEGER :: LC
        INTEGER :: i,j,k

        m    = 0
        UN   = 0
        TE   = 0
        HU   = 0
        TH   = 0
        TT   = 0
        LC   = 0

        fileoutput = 'pres_yzprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(17:17))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(16:16))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(15:15))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(14:14))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(13:13))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(12:12))

        OPEN(UNIT = 65, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        i = (NXL+NXR)/2
        
        do k = NZB,NZF
        do j = NYB,NYT

        write(65,*)p_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 65)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_ENERGY_YZ_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE energy,              ONLY : e_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 24) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
        REAL (KIND=8) :: x,y,z

        INTEGER :: m
        INTEGER :: UN
        INTEGER :: TE
        INTEGER :: HU
        INTEGER :: TH
        INTEGER :: TT
        INTEGER :: LC
        INTEGER :: i,j,k

        m    = 0
        UN   = 0
        TE   = 0
        HU   = 0
        TH   = 0
        TT   = 0
        LC   = 0

        fileoutput = 'ener_yzprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(17:17))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(16:16))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(15:15))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(14:14))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(13:13))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(12:12))

        OPEN(UNIT = 66, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        i = (NXL+NXR)/2
        
        do k = NZB,NZF
        do j = NYB,NYT

        write(66,*)e_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 66)

        END  SUBROUTINE
        !______________________________________________________________________        






       !______________________________________________________________________
        SUBROUTINE WRITE_VORTICITY_YZ_PROFILE(l)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE vorticity,           ONLY : omegax



        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 24) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
        REAL (KIND=8) :: x,y,z

        INTEGER :: m
        INTEGER :: UN
        INTEGER :: TE
        INTEGER :: HU
        INTEGER :: TH
        INTEGER :: TT
        INTEGER :: LC
        INTEGER :: i,j,k

        m    = 0
        UN   = 0
        TE   = 0
        HU   = 0
        TH   = 0
        TT   = 0
        LC   = 0

        fileoutput = 'vort_yzprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(17:17))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(16:16))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(15:15))
        m = (m-HU)/10
        TH = MOD(m,10)
        CALL VALUE(TH,fileoutput(14:14))
        m = (m-TH)/10
        TT = MOD(m,10)
        CALL VALUE(TT,fileoutput(13:13))
        m = (m-TT)/10
        LC = MOD(m,10)
        CALL VALUE(LC,fileoutput(12:12))

        OPEN(UNIT = 67, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        i = (NXL+NXR)/2

        do k = NZB,NZF
        do j = NYB,NYT

        write(67,*)omegax(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 67)

        END  SUBROUTINE
        !______________________________________________________________________





