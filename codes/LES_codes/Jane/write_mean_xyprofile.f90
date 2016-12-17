        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 02-02-2010        Last Modified on 17-05-2010

        !_____________________________________________________________________






        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_DENSITY_XY_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE continuity,          ONLY : rho_tot,rho_ts   !rho_mean !,rho_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 29) :: fileoutput
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

        fileoutput = 'den_mean_xyprof_000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(22:22))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(21:21))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(20:20))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(19:19))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(18:18))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(17:17))

        OPEN(UNIT = 91, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2

        do j = NYB,NYT
        do i = NXL,NXR
        !do j = NYB,NYT

        write(91,*)rho_tot(i,j,k), sqrt(rho_ts(i,j,k)) !rho_mean(i,j,k) !,rho_ms(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 91)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_U_VELOCITY_XY_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE x_momentum,          ONLY : u_tot,u_ts !u_mean !,u_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 29) :: fileoutput
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

        fileoutput = 'uvel_mean_xyprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(22:22))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(21:21))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(20:20))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(19:19))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(18:18))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(17:17))

        OPEN(UNIT = 92, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        k = (NZB+NZF)/2
        
        do j = NYB,NYT
        do i = NXL,NXR
        !do j = NYB,NYT

        write(92,*)u_tot(i,j,k), sqrt(u_ts(i,j,k)) !u_mean(i,j,k) !,u_ms(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 92)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_V_VELOCITY_XY_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE y_momentum,          ONLY : v_tot,v_ts  !v_mean !,v_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 29) :: fileoutput
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

        fileoutput = 'vvel_mean_xyprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(22:22))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(21:21))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(20:20))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(19:19))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(18:18))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(17:17))

        OPEN(UNIT = 93, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        k = (NZB+NZF)/2
        
        do j = NYB,NYT
        do i = NXL,NXR
        !do j = NYB,NYT

        write(93,*)v_tot(i,j,k),sqrt(v_ts(i,j,k))  !,v_mean(i,j,k) !,v_ms(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 93)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_W_VELOCITY_XY_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE z_momentum,          ONLY : w_tot,w_ts !w_mean !,w_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 29) :: fileoutput
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

        fileoutput = 'wvel_mean_xyprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(22:22))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(21:21))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(20:20))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(19:19))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(18:18))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(17:17))

        OPEN(UNIT = 94, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)
        
        k = (NZB+NZF)/2
        
        do i = NXL,NXR
        do j = NYB,NYT

        write(94,*) w_tot(i,j,k),sqrt(w_ts(i,j,k)) !w_mean(i,j,k) !,w_ms(i,j,k)

        end do
        end do

        CLOSE ( UNIT = 94)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_PRESSURE_XY_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE pressure,            ONLY : p_tot,p_ts ! p_mean,p_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 29) :: fileoutput
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

        fileoutput = 'pres_mean_xyprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(22:22))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(21:21))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(20:20))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(19:19))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(18:18))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(17:17))

        OPEN(UNIT = 95, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2
        
        do i = NXL,NXR
        do j = NYB,NYT
        write(95,*)p_tot(i,j,k),sqrt(p_ts(i,j,k)) !p_mean(i,j,k),p_ms(i,j,k)
        end do
        end do

        CLOSE ( UNIT = 95)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_MEAN_ENERGY_XY_PROFILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE energy,              ONLY : e_tot !e_mean !,e_ms

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 29) :: fileoutput
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

        fileoutput = 'ener_mean_xyprof000000.dat'
        m = l
        UN = MOD(m,10)
        CALL VALUE(UN,fileoutput(22:22))
        m = (m - UN)/10
        TE = MOD(m,10)
        CALL VALUE(TE,fileoutput(21:21))
        m = (m-TE)/10
        HU = MOD(m,10)
        CALL VALUE(HU,fileoutput(20:20))
        m = (m-HU)/10
        TH = MOD(m,10) 
        CALL VALUE(TH,fileoutput(19:19))
        m = (m-TH)/10
        TT = MOD(m,10) 
        CALL VALUE(TT,fileoutput(18:18))
        m = (m-TT)/10
        LC = MOD(m,10) 
        CALL VALUE(LC,fileoutput(17:17))

        OPEN(UNIT = 96, FILE = fileoutput , ACTION ='WRITE')
        120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        k = (NZB+NZF)/2
        
        do i = NXL,NXR
        do j = NYB,NYT
        write(96,*)e_tot(i,j,k) !e_mean(i,j,k) !,e_ms(i,j,k)
        end do
        end do

        CLOSE ( UNIT = 96)

        END  SUBROUTINE
        !______________________________________________________________________        







