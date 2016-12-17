        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 02-02-2010        Last Modified on 17-05-2010

        !_____________________________________________________________________






        !______________________________________________________________________
        SUBROUTINE WRITE_DENSITY_XZ_PROFILE(l)
       
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

        fileoutput = 'den_xzprof_000000.dat'
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

        OPEN(UNIT = 51, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        j = (NYB+NYT)/2

        do k = NZB,NZF
        do i = NXL,NXR

        write(51,*)rho_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 51)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_U_VELOCITY_XZ_PROFILE(l)
       
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

        fileoutput = 'uvel_xzprof000000.dat'
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

        OPEN(UNIT = 52, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        j = (NYB+NYT)/2
        
        do k = NZB,NZF
        do i = NXL,NXR

        write(52,*)u_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 52)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_V_VELOCITY_XZ_PROFILE(l)
       
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

        fileoutput = 'vvel_xzprof000000.dat'
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

        OPEN(UNIT = 53, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        j = (NYB+NYT)/2
        
        do k = NZB,NZF
        do i = NXL,NXR

        write(53,*)v_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 53)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_W_VELOCITY_XZ_PROFILE(l)
       
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

        fileoutput = 'wvel_xzprof000000.dat'
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

        OPEN(UNIT = 54, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        j = (NYB+NYT)/2
        
        do k = NZB,NZF
        do i = NXL,NXR

        write(54,*)w_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 54)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_PRESSURE_XZ_PROFILE(l)
       
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

        fileoutput = 'pres_xzprof000000.dat'
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

        OPEN(UNIT = 55, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        j = (NYB+NYT)/2
        
        do k = NZB,NZF
        do i = NXL,NXR

        write(55,*)p_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 55)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_ENERGY_XZ_PROFILE(l)
       
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

        fileoutput = 'ener_xzprof000000.dat'
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

        OPEN(UNIT = 56, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        j = (NYB+NYT)/2
        
        do k = NZB,NZF
        do i = NXL,NXR

        write(56,*)e_cur(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 56)

        END  SUBROUTINE
        !______________________________________________________________________        





        !______________________________________________________________________
        SUBROUTINE WRITE_VORTICITY_XZ_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE vorticity,           ONLY : omegay

       

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

        fileoutput = 'vort_xzprof000000.dat'
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

        OPEN(UNIT = 57, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        j = (NYB+NYT)/2
        
        do k = NZB,NZF
        do i = NXL,NXR

        write(57,*) omegay(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 57)

        END  SUBROUTINE
        !______________________________________________________________________        





