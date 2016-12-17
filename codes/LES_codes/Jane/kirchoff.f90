
        !_____________________________________________________________________
        
        ! PROGRAM WRIITEN BY SUBRAMANIAN GANESH

        ! Written on 05-07-2011        Last Modified on 05-07-2011

        !_____________________________________________________________________




        !______________________________________________________________________
        SUBROUTINE WRITE_KIRCHOFF(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF

        USE kirchoff        

       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 30) :: fileoutput
        REAL (KIND=8) :: x,y,z
        INTEGER :: i,j,k


        WRITE(fileoutput,'(a12,I6.6,a)') 'kirchoff_px_',l,'.dat'

        OPEN(301,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(301) p_px,dpdn_px
        CLOSE ( UNIT = 301)


        WRITE(fileoutput,'(a12,I6.6,a)') 'kirchoff_nx_',l,'.dat'

        OPEN(302,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(302) p_nx,dpdn_nx
        CLOSE ( UNIT = 302)




        WRITE(fileoutput,'(a12,I6.6,a)') 'kirchoff_py_',l,'.dat'

        OPEN(303,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(303) p_py,dpdn_py
        CLOSE ( UNIT = 303)


        WRITE(fileoutput,'(a12,I6.6,a)') 'kirchoff_ny_',l,'.dat'

        OPEN(304,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(304) p_ny,dpdn_ny
        CLOSE ( UNIT = 304)




        WRITE(fileoutput,'(a12,I6.6,a)') 'kirchoff_pz_',l,'.dat'

        OPEN(305,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(305) p_pz,dpdn_pz
        CLOSE ( UNIT = 305)


        WRITE(fileoutput,'(a12,I6.6,a)') 'kirchoff_nz_',l,'.dat'

        OPEN(306,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(306) p_nz,dpdn_nz
        CLOSE ( UNIT = 306)

        END  SUBROUTINE
        !______________________________________________________________________        





        !______________________________________________________________________
        SUBROUTINE CALCULATE_KIRCHOFF1

        USE GRID_DIMENSIONS

        USE pressure,        ONLY : p_x,p_y,p_z
        USE kirchoff

        IMPLICIT NONE

        INTEGER :: i,j,k

        do k = NZB,NZF
        do j = NYB,NYT
        dpdn_px(j,k) = p_x(NX_PXK,j,k)
        dpdn_nx(j,k) = 0.0d0 - p_x(NX_NXK,j,k)
        end do
        end do



        do k = NZB,NZF
        do i = NXL,NXR
        dpdn_py(i,k) = p_y(i,NY_PYK,k)
        dpdn_ny(i,k) = 0.0d0 - p_y(i,NY_NYK,k)
        end do
        end do


        do j = NYB,NYT
        do i = NXL,NXR
        dpdn_pz(i,j) = p_z(i,j,NZ_PZK)
        dpdn_nz(i,j) = 0.0d0 - p_z(i,j,NZ_NZK)
        end do
        end do


        END SUBROUTINE
        !______________________________________________________________________



        !______________________________________________________________________
        SUBROUTINE CALCULATE_KIRCHOFF2

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS, ONLY : NXB,NYSB,NYST,NZSB,NZSF
        USE pressure,        ONLY : p_cur,p_x,p_y,p_z
        USE kirchoff

        IMPLICIT NONE

        INTEGER :: i,j,k


        do k = NZB,NZF
        do j = NYB,NYT
        p_px(j,k) = p_cur(NX_PXK,j,k)
        p_nx(j,k) = p_cur(NX_NXK,j,k)

        dpdn_px(j,k) = 0.5d0 * (dpdn_px(j,k) + p_x(NX_PXK,j,k))
        dpdn_nx(j,k) = 0.5d0 * (dpdn_nx(j,k) - p_x(NX_NXK,j,k))
        end do
        end do



        do k = NZB,NZF
        do i = NXL,NXR
        p_py(i,k) = p_cur(i,NY_PYK,k)
        p_ny(i,k) = p_cur(i,NY_NYK,k)

        dpdn_py(i,k) = 0.5d0 * (dpdn_py(i,k) + p_y(i,NY_PYK,k))
        dpdn_ny(i,k) = 0.5d0 * (dpdn_ny(i,k) - p_y(i,NY_NYK,k))
        end do
        end do


        do j = NYB,NYT
        do i = NXL,NXR
        p_pz(i,j) = p_cur(i,j,NZ_PZK)
        p_nz(i,j) = p_cur(i,j,NZ_NZK)

        dpdn_pz(i,j) = 0.5d0 * (dpdn_pz(i,j) + p_z(i,j,NZ_PZK))
        dpdn_nz(i,j) = 0.5d0 * (dpdn_nz(i,j) - p_z(i,j,NZ_NZK))
        end do
        end do


        END SUBROUTINE
        !_______________________________________________________________________





        !______________________________________________________________________
        SUBROUTINE CALCULATE_MONITORPOINTS(l)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS, ONLY : NXB,NYSB,NYST,NZSB,NZSF
        USE continuity,      ONLY : rho_cur
        USE x_momentum,      ONLY : u_cur
        USE y_momentum,      ONLY : v_cur
        USE z_momentum,      ONLY : w_cur
        USE pressure,        ONLY : p_cur
        USE Monitor_Points
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: l
        INTEGER :: i,j,k,rem,fac

        rem = MOD(l,10)
        if (rem == 0 ) then
        fac = l/10 
        !write(*,*) fac

        do i = 1,NP_MP
        rho_mp(i,fac) = rho_cur(NX_MP(i),NY_MP(i),NZ_MP(i))
        !write(*,*) i,fac,rho_mp(i,fac)
        u_mp(i,fac)   = u_cur(NX_MP(i),NY_MP(i),NZ_MP(i))
        v_mp(i,fac)   = v_cur(NX_MP(i),NY_MP(i),NZ_MP(i))
        w_mp(i,fac)   = w_cur(NX_MP(i),NY_MP(i),NZ_MP(i))
        p_mp(i,fac)   = p_cur(NX_MP(i),NY_MP(i),NZ_MP(i))
        end do

        end if
        
        

        END SUBROUTINE
        !_______________________________________________________________________





        !______________________________________________________________________
        SUBROUTINE WRITE_MONITORPOINTS(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE Monitor_Points
       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 30) :: fileoutput
        REAL (KIND=8) :: x,y,z
        INTEGER :: i,j,k


        WRITE(fileoutput,'(a7,I6.6,a)') 'rho_mp_',l,'.dat'

        OPEN(301,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(301) rho_mp
        CLOSE ( UNIT = 301)


        WRITE(fileoutput,'(a5,I6.6,a)') 'u_mp_',l,'.dat'
        OPEN(302,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(302) u_mp
        CLOSE ( UNIT = 302)


        WRITE(fileoutput,'(a5,I6.6,a)') 'v_mp_',l,'.dat'
        OPEN(303,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(303) v_mp
        CLOSE ( UNIT = 303)


        WRITE(fileoutput,'(a5,I6.6,a)') 'w_mp_',l,'.dat'
        OPEN(304,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(304) u_mp
        CLOSE ( UNIT = 304)


        WRITE(fileoutput,'(a5,I6.6,a)') 'p_mp_',l,'.dat'
        OPEN(305,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(305) p_mp
        CLOSE ( UNIT = 305)


        END  SUBROUTINE
        !______________________________________________________________________        





        !______________________________________________________________________
        SUBROUTINE CALCULATE_MONITORLINES(l)

        USE GRID_DIMENSIONS, ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE GRID_DIMENSIONS, ONLY : NXB,NYSB,NYST,NZSB,NZSF
        USE continuity,      ONLY : rho_cur
        USE x_momentum,      ONLY : u_cur
        USE y_momentum,      ONLY : v_cur
        USE z_momentum,      ONLY : w_cur
        USE pressure,        ONLY : p_cur
        USE lighthill_source,ONLY : lh_s
        USE Monitor_Lines
        IMPLICIT NONE
        INTEGER, INTENT(INOUT) :: l
        INTEGER :: i,j,k,rem,fac

        rem = MOD(l,1)
        if (rem == 0 ) then
        fac = l/1
        !write(*,*) fac

        do j = 1,NL_ML
        do i = NXL,NXR
        rho_ml(j,i,fac)  = rho_cur(i,NY_ML(j),NZ_ML(j))
        !write(*,*) i,fac,rho_mp(i,fac)
        u_ml(j,i,fac)    = u_cur(i,NY_ML(j),NZ_ML(j))
        v_ml(j,i,fac)    = v_cur(i,NY_ML(j),NZ_ML(j))
        w_ml(j,i,fac)    = w_cur(i,NY_ML(j),NZ_ML(j))
        p_ml(j,i,fac)    = p_cur(i,NY_ML(j),NZ_ML(j))
        lh_s_ml(j,i,fac) = lh_s(i,NY_ML(j),NZ_ML(j))
        end do
        end do

        end if
        
        

        END SUBROUTINE
        !_______________________________________________________________________





        !______________________________________________________________________
        SUBROUTINE WRITE_MONITORLINES(l)
       
        USE GRID_DIMENSIONS,     ONLY : NXL,NXR,NYB,NYT,NZB,NZF
        USE Monitor_Lines
       

        IMPLICIT NONE
        INTEGER, INTENT(IN)          :: l

        CHARACTER(len = 40) :: fileoutput
        REAL (KIND=8) :: x,y,z
        INTEGER :: i,j,k


        WRITE(fileoutput,'(a7,I6.6,a)') 'rho_ml_',l,'.dat'

        OPEN(301,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(301) rho_ml
        CLOSE ( UNIT = 301)


        WRITE(fileoutput,'(a5,I6.6,a)') 'u_ml_',l,'.dat'
        OPEN(302,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(302) u_ml
        CLOSE ( UNIT = 302)


        WRITE(fileoutput,'(a5,I6.6,a)') 'v_ml_',l,'.dat'
        OPEN(303,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(303) v_ml
        CLOSE ( UNIT = 303)


        WRITE(fileoutput,'(a5,I6.6,a)') 'w_ml_',l,'.dat'
        OPEN(304,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(304) u_ml
        CLOSE ( UNIT = 304)


        WRITE(fileoutput,'(a5,I6.6,a)') 'p_ml_',l,'.dat'
        OPEN(305,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(305) p_ml
        CLOSE ( UNIT = 305)

        WRITE(fileoutput,'(a8,I6.6,a)') 'lh_s_ml_',l,'.dat'
        OPEN(306,FILE=fileoutput,form='unformatted') !,recl=NGP)
        write(306) lh_s_ml
        CLOSE ( UNIT = 306)


        END  SUBROUTINE
        !______________________________________________________________________        







