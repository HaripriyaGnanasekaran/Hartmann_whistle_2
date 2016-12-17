        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 09-11-2009        Last Modified on 10-05-2010

        !_____________________________________________________________________













        !______________________________________________________________________
        SUBROUTINE READ_FILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity,      ONLY : rho_cur
        USE x_momentum,      ONLY : u_cur
        USE y_momentum,      ONLY : v_cur
        USE z_momentum,      ONLY : w_cur
        USE pressure,        ONLY : p_cur
        USE energy,          ONLY : e_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
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

        fileoutput = 'acoustic_000000.dat'
        !fileoutput = 'bacustic_000000.dat'
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

        !OPEN(UNIT = 123, FILE = fileoutput , ACTION ='READ')
        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       
        OPEN(UNIT=123,FILE=fileoutput,form='unformatted')
 

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !read(123,*)u_cur(i,j,k),v_cur(i,j,k),w_cur(i,j,k),rho_cur(i,j,k),p_cur(i,j,k),e_cur(i,j,k)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 123)

        END  SUBROUTINE
        !______________________________________________________________________        










        !______________________________________________________________________
        SUBROUTINE READ_MEAN(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity,      ONLY : rho_tot !rho_mean
        USE x_momentum,      ONLY : u_tot   !u_mean
        USE y_momentum,      ONLY : v_tot   !v_mean
        USE z_momentum,      ONLY : w_tot   !w_mean
        USE pressure,        ONLY : p_tot   !p_mean
        USE energy,          ONLY : e_tot   !e_mean
        !USE vorticity,       ONLY : omega_mean

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
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

        fileoutput = 'meandata_000000.dat'
        !fileoutput = 'bacustic_000000.dat'
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

        OPEN(UNIT = 24, FILE = fileoutput , ACTION ='READ')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       

        do i = NXL,NXR
        do j = NYB,NYT
        do k = NZB,NZF

        !read(24,*)u_mean(i,j,k),v_mean(i,j,k),w_mean(i,j,k),rho_mean(i,j,k),p_mean(i,j,k),e_mean(i,j,k)!,omega_mean(i,j,k)

        end do
        end do
        end do

        CLOSE ( UNIT = 24)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE READ_MEAN_SQUARE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        !USE continuity,      ONLY : rho_ms
        !USE x_momentum,      ONLY : u_ms
        !USE y_momentum,      ONLY : v_ms
        !USE z_momentum,      ONLY : w_ms
        !USE pressure,        ONLY : p_ms
        !USE energy,          ONLY : e_ms
        !USE vorticity,       ONLY : omega_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
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

        fileoutput = 'mean_sqr_000000.dat'
        !fileoutput = 'bacustic_000000.dat'
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

        OPEN(UNIT = 25, FILE = fileoutput , ACTION ='READ')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       

        do i = NXL,NXR
        do j = NYB,NYT
        do k = NZB,NZF

        !read(25,*)p_ms(i,j,k)!,omega_ms(i,j,k)

        end do
        end do
        end do

        CLOSE ( UNIT = 25)

        END  SUBROUTINE
        !______________________________________________________________________        












        !______________________________________________________________________
        SUBROUTINE READ_DENSITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity,      ONLY : rho_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
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

        fileoutput = 'density__000000.dat'
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

        !OPEN(UNIT = 131, FILE = fileoutput , ACTION ='READ')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       
        OPEN(UNIT=131,FILE=fileoutput,form='unformatted')

        read(131) rho_cur
        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !read(131) (rho_cur(i,j,k),k=NZB,NZF)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 131)

        END  SUBROUTINE
        !______________________________________________________________________        













        !______________________________________________________________________
        SUBROUTINE READ_U_VELOCITY(l)
     
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE x_momentum,   ONLY : u_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
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

        fileoutput = 'uvelocity000000.dat'
        !fileoutput = 'bacustic_000000.dat'
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

        !OPEN(UNIT = 132, FILE = fileoutput , ACTION ='READ')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       
        OPEN(UNIT=132,FILE=fileoutput,form='unformatted')
        read(132) u_cur
        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !read(132) (u_cur(i,j,k),k=NZB,NZF)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 132)

        END  SUBROUTINE
        !______________________________________________________________________        













        !______________________________________________________________________
        SUBROUTINE READ_V_VELOCITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE y_momentum,   ONLY : v_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
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

        fileoutput = 'vvelocity000000.dat'
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

        !OPEN(UNIT = 133, FILE = fileoutput , ACTION ='READ')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       
        OPEN(UNIT=133,FILE=fileoutput,form='unformatted')
        read(133) v_cur

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !read(133) (v_cur(i,j,k),k=NZB,NZF)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 133)

        END  SUBROUTINE
        !______________________________________________________________________        













        !______________________________________________________________________
        SUBROUTINE READ_W_VELOCITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE z_momentum,   ONLY : w_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
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

        fileoutput = 'wvelocity000000.dat'
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

        !OPEN(UNIT = 134, FILE = fileoutput , ACTION ='READ')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       
        OPEN(UNIT=134,FILE=fileoutput,form='unformatted')
        read(134) w_cur

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !read(134) (w_cur(i,j,k),k=NZB,NZF)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 134)

        END  SUBROUTINE
        !______________________________________________________________________        













        !______________________________________________________________________
        SUBROUTINE READ_PRESSURE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE pressure,     ONLY : p_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
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

        fileoutput = 'pressure_000000.dat'
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

        !OPEN(UNIT = 135, FILE = fileoutput , ACTION ='READ')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       
        OPEN(UNIT=135,FILE=fileoutput,form='unformatted')
        read(135) p_cur

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !read(135) (p_cur(i,j,k),k=NZB,NZF)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 135)

        END  SUBROUTINE
        !______________________________________________________________________        













        !______________________________________________________________________
        SUBROUTINE READ_ENERGY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE energy,       ONLY : e_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
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

        fileoutput = 'energy_T_000000.dat'
        !fileoutput = 'bacustic_000000.dat'
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

        !OPEN(UNIT = 136, FILE = fileoutput , ACTION ='READ')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       
        OPEN(UNIT=136,FILE=fileoutput,form='unformatted')
        read(136) e_cur

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !read(136) (e_cur(i,j,k),k=NZB,NZF)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 136)

        END  SUBROUTINE
        !______________________________________________________________________        







        !______________________________________________________________________
        SUBROUTINE READ_VORTICITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE vorticity,       ONLY : omegax,omegay,omegaz

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
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

        fileoutput = 'vorticity000000.dat'
        !fileoutput = 'bacustic_000000.dat'
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

        !OPEN(UNIT = 137, FILE = fileoutput , ACTION ='READ')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
       
        OPEN(UNIT=137,FILE=fileoutput,form='unformatted')
        read(137) omegax,omegay,omegaz

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !read(137)(omegax(i,j,k),omegay(i,j,k),omegaz(i,j,k),k=NZB,NZF)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 137)

        END  SUBROUTINE
        !______________________________________________________________________        






        !______________________________________________________________________
        SUBROUTINE READ_LIGHTHILL_SOURCE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE lighthill_source

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
        INTEGER :: i,j,k


        fileoutput = 'lh_s_000000.dat'
        WRITE(fileoutput,'(a5,I6.6,a)') fileoutput,l,'.dat'

        !OPEN(UNIT = 35, FILE = fileoutput , ACTION ='READ')
        OPEN(UNIT = 135, FILE = fileoutput , form = 'unformatted' ) !ACTION ='READ')
        read(135) lh_s


        CLOSE ( UNIT = 135)

        END  SUBROUTINE
        !______________________________________________________________________        




