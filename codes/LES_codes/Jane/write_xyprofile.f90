        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 02-02-2010        Last Modified on 17-05-2010

        !_____________________________________________________________________






        !______________________________________________________________________
        SUBROUTINE WRITE_DENSITY_XY_PROFILE(l)
       
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

        fileoutput = 'den_xyprof_000000.dat'
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

        OPEN(UNIT = 41, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2

        do j = NYB,NYT
        do i = NXL,NXR

        write(41,*)rho_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 41)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_U_VELOCITY_XY_PROFILE(l)
       
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

        fileoutput = 'uvel_xyprof000000.dat'
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

        OPEN(UNIT = 42, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        k = (NZB+NZF)/2

        do j = NYB,NYT        
        do i = NXL,NXR

        write(42,*)u_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 42)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_V_VELOCITY_XY_PROFILE(l)
       
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

        fileoutput = 'vvel_xyprof000000.dat'
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

        OPEN(UNIT = 43, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        k = (NZB+NZF)/2

        do j = NYB,NYT        
        do i = NXL,NXR

        write(43,*)v_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 43)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_W_VELOCITY_XY_PROFILE(l)
       
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

        fileoutput = 'wvel_xyprof000000.dat'
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

        OPEN(UNIT = 44, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        k = (NZB+NZF)/2

        do j = NYB,NYT        
        do i = NXL,NXR

        write(44,*)w_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 44)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_PRESSURE_XY_PROFILE(l)
       
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

        fileoutput = 'pres_xyprof000000.dat'
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

        OPEN(UNIT = 45, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2

        do j = NYB,NYT        
        do i = NXL,NXR
        write(45,*)p_cur(i,j,k)
        end do
        end do

        CLOSE ( UNIT = 45)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_ENERGY_XY_PROFILE(l)
       
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

        fileoutput = 'ener_xyprof000000.dat'
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

        OPEN(UNIT = 46, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2
       
        do j = NYB,NYT 
        do i = NXL,NXR
        write(46,*)e_cur(i,j,k)
        end do
        end do

        CLOSE ( UNIT = 46)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_LIGHTHILL_SOURCE_XY_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE lighthill_source

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 24) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
        REAL (KIND=8) :: x
        REAL (KIND=8) :: y    
        INTEGER :: i,j,k


        fileoutput = 'lh_s_xyprof000000.dat'
        WRITE(fileoutput,'(a11,I6.6,a)') fileoutput,l,'.dat'

        OPEN(UNIT = 45, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2

        do j = NYB,NYT        
        do i = NXL,NXR
        write(45,*)lh_s(i,j,k)
        end do
        end do

        CLOSE ( UNIT = 45)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_TEMPERATURE_XY_PROFILE(l)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE energy



        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 24) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
        REAL (KIND=8) :: x
        REAL (KIND=8) :: y
        INTEGER :: i,j,k


        fileoutput = 'T_xyprof000000.dat'
        WRITE(fileoutput,'(a8,I6.6,a)') fileoutput,l,'.dat'

        OPEN(UNIT = 45, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2

        do j = NYB,NYT
        do i = NXL,NXR
        write(45,*)T(i,j,k)
        end do
        end do

        CLOSE ( UNIT = 45)

        END  SUBROUTINE
        !______________________________________________________________________        





        !______________________________________________________________________
        SUBROUTINE WRITE_VORTICITY_XY_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE vorticity,           ONLY : omegaz

       

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

        fileoutput = 'vort_xyprof000000.dat'
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

        OPEN(UNIT = 47, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2
       
        do j = NYB,NYT 
        do i = NXL,NXR
        write(47,*)omegaz(i,j,k)
        end do
        end do

        CLOSE ( UNIT = 47)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_SPL_XY_PROFILE(l)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE sound_field,         ONLY : spl_f



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

        fileoutput = 'spl__xyprof000000.dat'
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

        OPEN(UNIT = 48, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2

        do j = NYB,NYT
        do i = NXL,NXR
        write(48,*)spl_f(i,j,k)
        end do
        end do

        CLOSE ( UNIT = 48)

        END  SUBROUTINE
        !______________________________________________________________________        












