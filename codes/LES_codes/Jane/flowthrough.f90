! PROGRAM TO COMPUTE THE FLOWTHROUGH TIME AND NO OF STEPS NEEDED
! CODED ON 30 JUNE 2015
! by Ramanathan

!***********************************************************************************
SUBROUTINE DEVELOPMENTMODULE(l,t)

USE GRID_DIMENSIONS,             ONLY:NXO,NXR
USE GRID_SIZE_TIME_STEP,         ONLY:dx1,axo,dt
USE FLOW_PARAMETERS,             ONLY:v_jet

REAL(KIND=8) :: flow_through_time,xdistance,no_of_steps
INTEGER :: one_ft_time_steps
INTEGER,INTENT(OUT) :: t
INTEGER,INTENT(IN) :: l

!x-distance - gives the domain length
CALL X_COORD(xdistance,NXR)
xdistance = xdistance!NXO * dx1 + dx1*axo*((axo**(NXR-NXO) -1.0d0)/(axo-1.0d0))

!flow through time is computed using a approximate average jet velocity
!instead of integrating the axial velocity across the domain
flow_through_time = (xdistance)/(0.5*v_jet)

!number of numerical steps in time required for one flow through (approx.)
no_of_steps = flow_through_time/dt

!writing the data to user... to cross check the info.
WRITE(*,*)'**********************************************************************'
write(*,*) 'x-distance of domain                                  =',xdistance
WRITE(*,*) 'Time required per flow through                        =',flow_through_time
WRITE(*,*) 'No of time steps required to compute one flow through =',no_of_steps

one_ft_time_steps = int(no_of_steps)+1
tot_ft_time_steps = l*one_ft_time_steps
t=tot_ft_time_steps

WRITE(*,*) 'No of time steps required for FLOW TO DEVELOP         =',t

END SUBROUTINE
!********************************************************************************

SUBROUTINE FLOWTHROUGHMODULE(l,t)

USE GRID_DIMENSIONS,             ONLY:NXO,NXR,NX_st
USE GRID_SIZE_TIME_STEP,         ONLY:dx1,axo,dt
USE FLOW_PARAMETERS,             ONLY:v_jet

REAL(KIND=8) :: flow_through_time,xdistance,no_of_steps
INTEGER :: one_ft_time_steps
INTEGER,INTENT(OUT) :: t
INTEGER,INTENT(IN) :: l

CALL X_COORD(xdistance,NXR)
xdistance = xdistance
flow_through_time = (xdistance)/(0.5*v_jet)
no_of_steps = flow_through_time/dt
one_ft_time_steps = int(no_of_steps)+1
tot_ft_time_steps = l*one_ft_time_steps
t=tot_ft_time_steps
WRITE(*,*) 'No of time steps required for mean collection         =',t
WRITE(*,*) 'Duration for the mean (jet velocity*samples/distance) =',(v_jet*(t*dt))/xdistance,'(no unit)'
WRITE(*,*) '**********************************************************************'

END SUBROUTINE
!*******************************************************************************
