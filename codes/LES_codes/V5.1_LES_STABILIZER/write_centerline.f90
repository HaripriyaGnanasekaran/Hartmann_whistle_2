        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 02-02-2010        Last Modified on 17-05-2010

        !_____________________________________________________________________






        !______________________________________________________________________
        SUBROUTINE WRITE_DENSITY_CENTER_LINE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NY1,NY2

        USE continuity,          ONLY : rho_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
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

        fileoutput = 'den_cen__000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(15:15))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(14:14))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(13:13))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(12:12))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(11:11))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(10:10))

        OPEN(UNIT = 51, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2
        j = (NY1+NY2)/2

        do i = NXL,NXR

        write(51,*)rho_cur(i,j,k)

        end do

        CLOSE ( UNIT = 51)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_PRESSURE_CENTER_LINE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NY1,NY2

        USE pressure,            ONLY : p_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
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

        fileoutput = 'pres_cen_000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(15:15))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(14:14))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(13:13))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(12:12))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(11:11))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(10:10))

        OPEN(UNIT = 52, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2
        j = (NY1+NY2)/2
        
        do i = NXL,NXR
        write(52,*)p_cur(i,j,k)
        end do

        CLOSE ( UNIT = 52)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_ENERGY_CENTER_LINE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY : NY1,NY2

        USE energy,              ONLY : e_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
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

        fileoutput = 'ener_cen_000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(15:15))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(14:14))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(13:13))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(12:12))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(11:11))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(10:10))

        OPEN(UNIT = 53, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2
        j = (NY1+NY2)/2
        
        do i = NXL,NXR
        write(53,*)e_cur(i,j,k)
        end do

        CLOSE ( UNIT = 53)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_U_VELOCITY_CENTER_LINE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE x_momentum,          ONLY : u_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
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

        fileoutput = 'uvel_cen_000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(15:15))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(14:14))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(13:13))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(12:12))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(11:11))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(10:10))

        OPEN(UNIT = 54, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        k = (NZB+NZF)/2
        j = (NYB+NYT)/2
        
        do i = NXL,NXR

        write(54,*)u_cur(i,j,k)

        end do

        CLOSE ( UNIT = 54)

        END  SUBROUTINE
        !______________________________________________________________________        







