        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 20-08-2009        Last Modified on 10-05-2010

        !_____________________________________________________________________


		subroutine gen_data_for_lcs(l)
		
		USE grid_dimensions
        USE GRID_SIZE_TIME_STEP
        USE continuity
        USE x_momentum
        USE y_momentum
		
		implicit none
		integer :: i,j,k=0
		integer, intent(in) :: l
		real(kind=8) :: x,y
		character(len = 50) :: fhandle
		
		write(fhandle,'(a8,I5.5,a4)' ) "velocity",l, ".txt"
		open(unit = 864, file=fhandle, form='formatted')
		
		do i=nxl,nxr
			call x_coord(x,i)
			do j=nyb,nyt
				call y_coord(y,j)
				write(864,'(F7.2,F7.2,F7.5,F7.5)') x, y, u_cur(i,j,k), v_cur(i,j,k)
			end do
		end do
		
		
		end subroutine gen_data_for_lcs










        !______________________________________________________________________
        SUBROUTINE WRITE_FILE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF !,NGP

        USE continuity,   ONLY : rho_cur
        USE x_momentum,   ONLY : u_cur
        USE y_momentum,   ONLY : v_cur
        USE z_momentum,   ONLY : w_cur
        USE pressure,     ONLY : p_cur
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
        !INTEGER :: NP

        m    = 0
        UN   = 0
        TE   = 0
        HU   = 0
        TH   = 0
        TT   = 0
        LC   = 0
        !NP   = (NXR-NXL+1)*(NYT-NYB+1)*(NZF-NZB+1)*48
        !write(*,*) NP,(NXR-NXL+1)*(NYT-NYB+1)*(NZF-NZB+1) 

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

        !OPEN(UNIT = 12, FILE = fileoutput , ACTION ='WRITE')

        !OPEN(UNIT=12,FILE=fileoutput,form='unformatted',access='direct',recl=(NP*6*8))
        OPEN(UNIT=12,FILE=fileoutput,form='unformatted')! ,recl=(NGP))
        write(12,rec=6) u_cur,v_cur,w_cur,rho_cur,p_cur,e_cur

        CLOSE (UNIT=12)

        write(*,*) rho_cur(NXL,NYB,NZB)

        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF

        !write(12,*) u_cur(i,j,k),v_cur(i,j,k),w_cur(i,j,k),rho_cur(i,j,k),p_cur(i,j,k),e_cur(i,j,k)

        !end do
        !end do
        !end do

        
        !CLOSE ( UNIT = 12)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_P(NP,u,v,rho,p,e,dt)
        IMPLICIT NONE
       
        INTEGER, INTENT(IN) :: NP
        REAL (KIND=8), DIMENSION(1:NP), INTENT(IN) :: rho
        REAL (KIND=8), DIMENSION(1:NP), INTENT(IN) :: u
        REAL (KIND=8), DIMENSION(1:NP), INTENT(IN) :: v
        REAL (KIND=8), DIMENSION(1:NP), INTENT(IN) :: p
        REAL (KIND=8), DIMENSION(1:NP), INTENT(IN) :: e
        REAL (KIND=8),                INTENT(IN) :: dt

        REAL (KIND=8), PARAMETER :: pi = 3.14159265359
        REAL (KIND=8) :: t
        INTEGER :: i

        OPEN(UNIT = 13, FILE = 'spectrum_008000.dat' , ACTION ='WRITE')

        i = 1
        do
        if (i > NP) exit
        t = (i-1) * dt * 10.0D0
        write(13,*)t,rho(i),u(i),v(i),p(i),e(i)
        write(*,*) t,rho(i),u(i),v(i),p(i),e(i)
        i = i + 1
        end do
        CLOSE ( UNIT = 13)

        END  SUBROUTINE
        !______________________________________________________________________















        !______________________________________________________________________
        SUBROUTINE WRITE_DENSITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF !,NGP
        USE GRID_SIZE_TIME_STEP, ONLY : dx1,dy1,dz1,dx,dy,dz

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

        OPEN(121,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(121) rho_cur

        !OPEN(UNIT = 121, FILE = fileoutput , ACTION ='WRITE')
        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(121) (rho_cur(i,j,k),k=NZB,NZF,j=NYB,NYT)
        !end do
        !end do
        !end do
        CLOSE ( UNIT = 121)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_U_VELOCITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF !,NGP

        USE x_momentum,   ONLY : u_cur

       

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

        fileoutput = 'uvelocity000000.dat'
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

        OPEN(122,FILE=fileoutput,form='unformatted')!,recl=NGP)
        write(122) u_cur

        !OPEN(UNIT = 122, FILE = fileoutput , ACTION ='WRITE')
        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(122) (u_cur(i,j,k),k=NZB,NZF,j=NYB,NYT)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 122)

        END  SUBROUTINE
        !______________________________________________________________________        










        !______________________________________________________________________
        SUBROUTINE WRITE_V_VELOCITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF !,NGP

        USE y_momentum,   ONLY : v_cur

       

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

        OPEN(123,FILE=fileoutput,form='unformatted')!,recl=NGP)
        write(123) v_cur

        !OPEN(UNIT = 123, FILE = fileoutput , ACTION ='WRITE')
        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(123) (v_cur(i,j,k),k=NZB,NZF,j=NYB,NYT)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 123)

        END  SUBROUTINE
        !______________________________________________________________________        









        !______________________________________________________________________
        SUBROUTINE WRITE_W_VELOCITY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF !,NGP

        USE z_momentum,   ONLY : w_cur

       

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

        OPEN(124,FILE=fileoutput,form='unformatted')!,recl=NGP)
        write(124) w_cur

        !OPEN(UNIT = 124, FILE = fileoutput , ACTION ='WRITE')
        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(124) (w_cur(i,j,k),k=NZB,NZF,j=NYB,NYT)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 124)

        END  SUBROUTINE
        !______________________________________________________________________        












        !______________________________________________________________________
        SUBROUTINE WRITE_PRESSURE(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF !,NGP

        USE pressure,     ONLY : p_cur

       

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

        OPEN(125,FILE=fileoutput,form='unformatted')!,recl=NGP)
        write(125) p_cur

        !OPEN(UNIT = 125, FILE = fileoutput , ACTION ='WRITE')
        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(125) (p_cur(i,j,k),k=NZB,NZF,j=NYB,NYT)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 125)

        END  SUBROUTINE
        !______________________________________________________________________        











        !______________________________________________________________________
        SUBROUTINE WRITE_ENERGY(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF !,NGP

        USE energy,       ONLY : e_cur

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
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

        fileoutput = 'energy_T_000000.dat'
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

        OPEN(126,FILE=fileoutput,form='unformatted')!,recl=NGP)
        write(126) e_cur

        !OPEN(UNIT = 126, FILE = fileoutput , ACTION ='WRITE')
        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(126) (e_cur(i,j,k),k=NZB,NZF,j=NYB,NYT)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 126)

        END  SUBROUTINE
        !______________________________________________________________________        








        !______________________________________________________________________
        SUBROUTINE WRITE_VORTICITY(l)

        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF !,NGP

        USE vorticity,       ONLY :omegax,omegay,omegaz



        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 22) :: fileoutput
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


        OPEN(UNIT=127,FILE=fileoutput,form='unformatted')!,recl=NGP)
        write(127) omegax,omegay,omegaz

        !OPEN(UNIT = 127, FILE = fileoutput , ACTION ='WRITE')
        !120 FORMAT(F15.12,T18,F15.12,T36,F15.12,T54,F15.12,T78,F15.12)!,T100,F15.9,T120,F15.9)

        !do i = NXL,NXR
        !do j = NYB,NYT
        !do k = NZB,NZF
        !write(127) (omegax(i,j,k),omegay(i,j,k),omegaz(i,j,k),k=NZB,NZF)
        !end do
        !end do
        !end do

        CLOSE ( UNIT = 127)

        END  SUBROUTINE
        !______________________________________________________________________








        !______________________________________________________________________
        SUBROUTINE WRITE_LIGHTHILL_SOURCE(l)
       
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
        OPEN(UNIT = 135, FILE = fileoutput , form = 'unformatted' ) !ACTION ='WRITE')
        write(135) lh_s


        CLOSE ( UNIT = 135)

        END  SUBROUTINE
        !______________________________________________________________________        




