        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 09-12-2009        Last Modified on 10-05-2010

        !_____________________________________________________________________















        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity,   ONLY : rho_tot !,rho_mean
        USE x_momentum,   ONLY : u_tot   !,u_mean
        USE y_momentum,   ONLY : v_tot   !,v_mean
        USE z_momentum,   ONLY : w_tot   !,w_mean
        USE pressure,     ONLY : p_tot   !,p_mean
        USE energy,       ONLY : e_tot   !,e_mean
        !USE vorticity,    ONLY : omega_mean

       

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

        fileoutput = 'meandata_005000.dat'
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

        OPEN(UNIT = 15, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        !120 FORMAT(2F15.12,3F15.12)

        do i = NXL,NXR
        do j = NYB,NYT
        do k = NZB,NZF

        !write(15,*)u_mean(i,j,k),v_mean(i,j,k),w_mean(i,j,k),rho_mean(i,j,k),p_mean(i,j,k),e_mean(i,j,k) !,omega_mean(i,j,k)

        end do
        end do
        end do

        CLOSE ( UNIT = 15)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_SQUARE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity,   ONLY : rho_ts
        USE x_momentum,   ONLY : u_ts
        USE y_momentum,   ONLY : v_ts
        USE z_momentum,   ONLY : w_ts
        USE pressure,     ONLY : p_ts
        !USE energy,       ONLY : e_ms
        !USE vorticity,    ONLY : omega_ms

       

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

        OPEN(UNIT = 16, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        !120 FORMAT(2F15.12,3F15.12)

        do i = NXL,NXR
        do j = NYB,NYT
        do k = NZB,NZF

        !write(16,*)p_ms(i,j,k)

        end do
        end do
        end do

        CLOSE ( UNIT = 16)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_SPL(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE sound_field,  ONLY : spl_f

       

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

        fileoutput = 'splfield_000000.dat'
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

        OPEN(UNIT = 117, FILE = fileoutput , form = 'unformatted')!,ACTION ='WRITE')
        !OPEN(UNIT = 117, FILE = fileoutput , ACTION ='WRITE')
        
        write(117) spl_f
        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF

        !write(17,*) spl_f(i,j,k)

        !end do
        !end do
        !end do

        CLOSE ( UNIT = 117)

        END  SUBROUTINE
        !______________________________________________________________________        













        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_DENSITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity,      ONLY : rho_tot,rho_ts  !rho_mean !,rho_ms

       

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

        fileoutput = 'mean_den_000000.dat'
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

        !OPEN(UNIT = 131, FILE = fileoutput , ACTION ='WRITE')
        OPEN(UNIT = 131, FILE = fileoutput , form = 'unformatted') ! ACTION ='WRITE')
        write(131) rho_tot,rho_ts

        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(31,*) rho_mean(i,j,k)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 131)

        END  SUBROUTINE
        !______________________________________________________________________        










        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_U_VELOCITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE x_momentum,   ONLY : u_tot,u_ts  ! u_mean !,u_ms

       

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

        fileoutput = 'meanuvel_000000.dat'
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

        !OPEN(UNIT = 132, FILE = fileoutput , ACTION ='WRITE')
        OPEN(UNIT = 132, FILE = fileoutput , form = 'unformatted' ) !ACTION ='WRITE')
        write(132) u_tot,u_ts

        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(32,*) u_mean(i,j,k) !,u_ms(i,j,k)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 132)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_V_VELOCITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE y_momentum,   ONLY : v_tot,v_ts ! v_mean !,v_ms

       

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

        fileoutput = 'meanvvel_000000.dat'
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

        !OPEN(UNIT = 133, FILE = fileoutput , ACTION ='WRITE')
        OPEN(UNIT = 133, FILE = fileoutput , form = 'unformatted' ) !ACTION ='WRITE')
        write(133) v_tot,v_ts

        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(33,*) v_mean(i,j,k) !,v_ms(i,j,k)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 133)

        END  SUBROUTINE
        !______________________________________________________________________        










        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_W_VELOCITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE z_momentum,   ONLY : w_tot,w_ts !,w_mean !,w_ms

       

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

        fileoutput = 'meanwvel_000000.dat'
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

        !OPEN(UNIT = 134, FILE = fileoutput , ACTION ='WRITE')
        OPEN(UNIT = 134, FILE = fileoutput , form = 'unformatted') !ACTION ='WRITE')
        write(134) w_tot,w_ts 

        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(34,*) w_mean(i,j,k) !,w_ms(i,j,k)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 134)

        END  SUBROUTINE
        !______________________________________________________________________        












        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_PRESSURE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE pressure,     ONLY : p_tot,p_ts  !p_mean,p_ms

       

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

        fileoutput = 'meanpres_000000.dat'
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

        !OPEN(UNIT = 35, FILE = fileoutput , ACTION ='WRITE')
        OPEN(UNIT = 135, FILE = fileoutput , form = 'unformatted' ) !ACTION ='WRITE')
        write(135) p_tot,p_ts  !p_mean,p_ms

        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(35,*) p_mean(i,j,k),p_ms(i,j,k)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 135)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_ENERGY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE energy,       ONLY : e_tot  !e_mean !,e_ms

       

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

        fileoutput = 'mean_e_T_000000.dat'
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

        !OPEN(UNIT = 136, FILE = fileoutput , ACTION ='WRITE')
        OPEN(UNIT = 136, FILE = fileoutput , form = 'unformatted') !ACTION ='WRITE')
        write(136) e_tot   ! e_mean

        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(36,*) e_mean(i,j,k) !,e_ms(i,j,k)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 136)

        END  SUBROUTINE
        !______________________________________________________________________        







        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_CROSS_MOMENTUM(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE cross_momentum , ONLY : uv_tot,uw_tot,vw_tot

       

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

        fileoutput = 'cross_mom000000.dat'
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

        !OPEN(UNIT = 134, FILE = fileoutput , ACTION ='WRITE')
        OPEN(UNIT = 134, FILE = fileoutput , form = 'unformatted') !ACTION ='WRITE')
        write(134) uv_tot,uw_tot,vw_tot 

        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(34,*) w_mean(i,j,k) !,w_ms(i,j,k)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 134)

        END  SUBROUTINE
        !______________________________________________________________________        













        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_VORTICITY(l)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        !USE vorticity,       ONLY :omega_mean,omega_ms



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

        fileoutput = 'vorticity000000.dat'
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


        OPEN(UNIT = 37, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        do i = NXL,NXR
        do j = NYB,NYT
        do k = NZB,NZF
        !write(37,*) omega_mean(i,j,k),omega_ms(i,j,k)
        end do
        end do
        end do

        CLOSE ( UNIT = 37)

        END  SUBROUTINE
        !______________________________________________________________________







        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_LIGHTHILL_SOURCE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE lighthill_source

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
        INTEGER :: i,j,k


        fileoutput = 'mean_lh_s_000000.dat'
        WRITE(fileoutput,'(a10,I6.6,a)') fileoutput,l,'.dat'

        !OPEN(UNIT = 35, FILE = fileoutput , ACTION ='READ')
        OPEN(UNIT = 135, FILE = fileoutput , form = 'unformatted' ) !ACTION ='WRITE')
        write(135) lh_s_tot,lh_s_ts


        CLOSE ( UNIT = 135)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_RHOU_MOMENTUM(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE x_momentum,   ONLY : rhou_tot  ! u_mean !,u_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 30) :: fileoutput

        INTEGER :: i,j,k


        fileoutput = 'mean_rhou_mom_000000.dat'

        WRITE(fileoutput,'(a14,I6.6,a)') fileoutput,l,'.dat'


        OPEN(UNIT = 142, FILE = fileoutput , form = 'unformatted' ) !ACTION ='WRITE')
        write(142) rhou_tot


        CLOSE ( UNIT = 142)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_RHOV_MOMENTUM(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE y_momentum,   ONLY : rhov_tot  ! u_mean !,u_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 30) :: fileoutput

        INTEGER :: i,j,k


        fileoutput = 'mean_rhov_mom_000000.dat'

        WRITE(fileoutput,'(a14,I6.6,a)') fileoutput,l,'.dat'


        OPEN(UNIT = 143, FILE = fileoutput , form = 'unformatted' ) !ACTION ='WRITE')
        write(143) rhov_tot


        CLOSE ( UNIT = 143)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_RHOW_MOMENTUM(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE z_momentum,   ONLY : rhow_tot  ! u_mean !,u_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 30) :: fileoutput

        INTEGER :: i,j,k


        fileoutput = 'mean_rhow_mom_000000.dat'

        WRITE(fileoutput,'(a14,I6.6,a)') fileoutput,l,'.dat'


        OPEN(UNIT = 144, FILE = fileoutput , form = 'unformatted' ) !ACTION ='WRITE')
        write(144) rhow_tot


        CLOSE ( UNIT = 144)

        END  SUBROUTINE
        !______________________________________________________________________        






        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_TURB_KE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE cross_momentum , ONLY : turb_ke_tot
       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 30) :: fileoutput

        INTEGER :: i,j,k


        fileoutput = 'mean_turb_ke_000000.dat'

        WRITE(fileoutput,'(a13,I6.6,a)') fileoutput,l,'.dat'


        OPEN(UNIT = 145, FILE = fileoutput , form = 'unformatted' ) !ACTION ='WRITE')
        write(145) turb_ke_tot


        CLOSE ( UNIT = 145)

        END  SUBROUTINE
        !______________________________________________________________________        



