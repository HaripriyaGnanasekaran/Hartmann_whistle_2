        ! This program contains the modules which contains the main
        ! parameters
        !______________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 13-07-2009        Last Modified on 04-05-2010

        !______________________________________________________________________



        !______________________________________________________________________
        MODULE GRID_DIMENSIONS
        IMPLICIT NONE
        SAVE
        
        INTEGER, PARAMETER :: NXL = 0 
        INTEGER, PARAMETER :: NXR = 114!NUMBER OF GRID POINTS ALONG X AXIS
        INTEGER, PARAMETER :: NYB = -55
        INTEGER, PARAMETER :: NYT = 55!NUMBER OF GRID POINTS ALONG Y AXIS
        INTEGER, PARAMETER :: NZB = -55
        INTEGER, PARAMETER :: NZF = 55 !NUMBER OF GRID POINTS ALONG Z AXIS

        INTEGER, PARAMETER :: NXO = 100
        INTEGER, PARAMETER :: NXB  = 100
        INTEGER, PARAMETER :: NX_st = 32

        INTEGER, PARAMETER :: NY1 = -8  !JET EXIT
        INTEGER, PARAMETER :: NY2 = 8  !JET EXIT

        INTEGER, PARAMETER :: NYSB = -24
        INTEGER, PARAMETER :: NYST = 24

        INTEGER, PARAMETER :: NZSB = -24
        INTEGER, PARAMETER :: NZSF = 24

        INTEGER, PARAMETER :: NX_NXK = 3
        INTEGER, PARAMETER :: NX_PXK = 80

        INTEGER, PARAMETER :: NY_NYK = -42
        INTEGER, PARAMETER :: NY_PYK = 42

        INTEGER, PARAMETER :: NZ_NZK = -42
        INTEGER, PARAMETER :: NZ_PZK = 42
    
        !INTEGER, PARAMETER :: NGP = (NXR-NXL+1)*(NYT-NYB+1)*(NZF-NZB+1)
        !Total Number of grid points * 8 = Record length

        END MODULE GRID_DIMENSIONS
        !______________________________________________________________________


        !______________________________________________________________________
        MODULE GRID_SIZE_TIME_STEP
        USE GRID_DIMENSIONS
        IMPLICIT NONE
        SAVE
        
        REAL (KIND=8), PARAMETER :: dx1=1.25*0.02402355516!1.354521203624431d0*10.0d0**(-3)!4.9053086d0*10.0d0**(-4)
        REAL (KIND=8), PARAMETER :: dy1= 0.02402355516 !4.9053086d0*10.0d0**(-4)
        REAL (KIND=8), PARAMETER :: dz1 = dy1 !4.9053086d0*10.0d0**(-4)

        REAL (KIND=8)  :: dx
        REAL (KIND=8)  :: dy
        REAL (KIND=8)  :: dz


        REAL (KIND=8), DIMENSION(NXL:NXR-1) :: dxi
        REAL (KIND=8), DIMENSION(NYB:NYT-1) :: dyi
        REAL (KIND=8), DIMENSION(NZB:NZF-1) :: dzi

        REAL (KIND=8), PARAMETER :: axs = 1.01d0

        REAL (KIND=8), PARAMETER :: axo = 1.05d0
        REAL (KIND=8), PARAMETER :: ayo = 1.02d0
        REAL (KIND=8), PARAMETER :: azo = 1.02d0

        REAL (KIND=8), PARAMETER :: dt =  0.075d0*dy1/340.0D0
        REAL (KIND=8), PARAMETER :: jet_d = 16.0d0*dy1
        REAL (KIND=8), PARAMETER :: jet_r = 8.0d0*dy1
        
        REAL (KIND=8), PARAMETER :: d_eta_x = dy1
        REAL (KIND=8), PARAMETER :: d_eta_x_inv = 1.0D0/d_eta_x
        REAL (KIND=8), PARAMETER :: d_eta_y = dy1,d_eta_y_inv=1.0d0/d_eta_y
        REAL (KIND=8), PARAMETER :: d_eta_z = dy1,d_eta_z_inv=1.0d0/d_eta_z
                
        END MODULE GRID_SIZE_TIME_STEP
        !______________________________________________________________________

        
        !______________________________________________________________________
        MODULE FLOW_PARAMETERS
        USE GRID_SIZE_TIME_STEP
        IMPLICIT NONE
        SAVE
        
        REAL (KIND=8) :: pi   = 3.14159265359

        REAL (KIND=8), PARAMETER :: rho_atm = 1.25D0
        REAL (KIND=8), PARAMETER :: p_atm = 101325.0D0

        REAL (KIND=8), PARAMETER :: rho_c = 0.0220678d0
        REAL (KIND=8), PARAMETER :: p_c = 1823.85d0	
        REAL (KIND=8), PARAMETER :: T_c = 287.93d0

        REAL (KIND=8), PARAMETER :: rho_in= 0.02564d0
        REAL (KIND=8), PARAMETER :: p_in = 1823.85d0     
        REAL (KIND=8), PARAMETER :: T_in = 247.78853d0


        REAL (KIND=8), PARAMETER :: gama = 1.4D0  
        REAL (KIND=8), PARAMETER :: R = 287.04d0
        REAL (KIND=8), PARAMETER :: Cv = 717.6d0
        REAL (KIND=8), PARAMETER :: Cp = 1004.64d0
      
        REAL (KIND=8) :: c = 315.5555538!315.5555 !sqrt(1.4d0*101325.0d0/1.25d0)
       
        REAL (KIND=8), PARAMETER :: v_jet = 50.90722522d0
        REAL (KIND=8) :: mu = 1.587697d0*10.0d0**(-5) !1.983*10.0D0**(-5)
        REAL (KIND=8), PARAMETER :: Pr = 0.72d0
        REAL (KIND=8) :: kT = 0.024968688 !0.02766946d0

        END MODULE FLOW_PARAMETERS
        !______________________________________________________________________
        


        !______________________________________________________________________
        MODULE HIXON_FORWARD_PARAMETERS
        IMPLICIT NONE
        SAVE

        REAL (KIND=8) :: a = 0.21132486540518708 !(1.0D0 - 1.0D0/sqrt(3.0D0))*0.5D0
        REAL (KIND=8) :: ad = 0.78867513459481287 !1.0D0 - (1.0D0 - 1.0D0/sqrt(3.0D0))*0.5D0
        REAL (KIND=8) :: adi = 1.0d0/0.78867513459481287d0

        REAL (KIND=8) :: b1 = 2.9012462146853033 !25.0D0/12.0D0+17.0D0/(12.0D0*sqrt(3.0D0))
        REAL (KIND=8) :: b2 = -6.40562612162344 ! 0.0D0 - (4.0D0 + 25.0D0/ (6.0D0 * sqrt(3.0D0)))
        REAL (KIND=8) :: b3 = 5.598076211353316 !3.0D0 + 3.0D0*sqrt(3.0D0)*0.5D0
        REAL (KIND=8) :: b4 = -2.5842589165775225 ! 0.0D0 - (4.0D0/3.0D0 + 13.0D0/(6.0D0 * sqrt(3.0D0)))
        REAL (KIND=8) :: b5 = 0.49056261216234409 !0.25D0 + 5.0D0/(12.0D0 * sqrt(3.0D0) )

        END MODULE HIXON_FORWARD_PARAMETERS
        !______________________________________________________________________




        !______________________________________________________________________
        MODULE HIXON_BACKWARD_PARAMETERS
        IMPLICIT NONE
        SAVE

        REAL (KIND=8) :: a  =  0.21132486540518708 ! (1.0D0 - 1.0D0/sqrt(3.0D0))*0.5D0
        REAL (KIND=8) :: ad =  0.78867513459481287 ! 1.0D0 - (1.0D0 - 1.0D0/sqrt(3.0D0))*0.5D0
        REAL (KIND=8) :: adi = 1.0d0/0.78867513459481287d0

        REAL (KIND=8) :: b1 = -2.9012462146853033 ! 0.0D0 - (25.0D0/12.0D0 + 17.0D0/(12.0D0 * sqrt(3.0D0)))
        REAL (KIND=8) :: b2 = 6.40562612162344 ! 4.0D0 + 25.0D0/ (6.0D0 * sqrt(3.0D0))
        REAL (KIND=8) :: b3 = -5.598076211353316 ! 0.0D0 - (3.0D0 + 3.0D0*sqrt(3.0D0)*0.5D0)
        REAL (KIND=8) :: b4 = 2.5842589165775225 ! 4.0D0/3.0D0 + 13.0D0/(6.0D0 * sqrt(3.0D0))
        REAL (KIND=8) :: b5  = -0.49056261216234409 !0.0D0 - (0.25D0 + 5.0D0/(12.0D0 * sqrt(3.0D0) ))

        END MODULE HIXON_BACKWARD_PARAMETERS
        !______________________________________________________________________

