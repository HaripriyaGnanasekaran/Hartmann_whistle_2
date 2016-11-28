        !______________________________________________________________________

        ! Program written by SUBRAMANIAN GANESH
        
        ! Written on 12-01-2010        Last Modified on 17-05-2010

        !______________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE BUFFER_FILTER(u_cur)

        USE GRID_DIMENSIONS,     ONLY: NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS,     ONLY: NXB

        USE GRID_SIZE_TIME_STEP, ONLY:dx1,dy1,dz1
        USE BUFFER_FILTER_PARAMETERS, ONLY: x_buf

        IMPLICIT NONE

        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), INTENT(INOUT) :: u_cur
        REAL (KIND=8), DIMENSION(NXL:NXR,NYB:NYT,NZB:NZF), SAVE  :: u_fil

        REAL (KIND=8), PARAMETER :: alfa = 0.2d0
        REAL (KIND=8), PARAMETER :: beta = 1.5d0
        REAL (KIND=8), PARAMETER :: d0 = 1.5d0
        REAL (KIND=8), PARAMETER :: d1 = 0.0d0 - 0.25d0
        REAL (KIND=8)            :: u

        INTEGER :: i,j,k
        
        CALL COPY(u_cur,u_fil)

    
        
        do k = NZB+1,NZF-1
        do j = NYB+1,NYT-1
        do i = NXB+1,NXR-1

        u = u_cur(i-1,j,k) + u_cur(i+1,j,k)
        u = u + u_cur(i,j-1,k) + u_cur(i,j+1,k)
        u = u + u_cur(i,j,k-1) + u_cur(i,j,k+1)
        u = u*d1
        u = u + d0 * u_cur(i,j,k)

        u_fil(i,j,k)  = u_cur(i,j,k) - alfa*u*x_buf(i)
 
        end do
        end do
        end do

        !k = NZB
        !
        !do j = NYB+1,NYT-1
        !do i = NXB+1,NXR-1
        !
        !u = u_cur(i-1,j,k) + u_cur(i+1,j,k)
        !u = u + u_cur(i,j-1,k) + u_cur(i,j+1,k)
        !u = u + u_cur(i,j,NZF-1) + u_cur(i,j,NZB+1)
        !u = u*d1
        !u = u + d0 * u_cur(i,j,k)
        !
        !
        !u_fil(i,j,k)  = u_cur(i,j,k) - alfa*u*x_buf(i)
        !
        !u_fil(i,j,NZF) = u_fil(i,j,NZB)
        !end do
        !end do



        CALL COPY(u_fil,u_cur)

        !write(*,*) 'Buffer Done'

        END SUBROUTINE
        !______________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE INITIALIZE_BUFFER_FILTER
        
        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NXB
        USE BUFFER_FILTER_PARAMETERS

        IMPLICIT NONE

        INTEGER :: i
        REAL(KIND=8) :: x1,xb,xbw
        REAL(KIND=8) :: beta = 1.5d0
        
        CALL X_COORD(xb,NXB)
        CALL X_COORD(xbw,NXR)
        xbw = xbw - xb

        do i = NXB+1,NXR-1
        CALL X_COORD(x1,i)
        x1 = x1 - xb

        x_buf(i) = (x1/xbw)**beta
        end do
        
        

        END SUBROUTINE
        !______________________________________________________________________






