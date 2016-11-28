        SUBROUTINE CHECK_NRBC(ts)

        USE GRID_DIMENSIONS
        USE GRID_SIZE_TIME_STEP
        USE FLOW_PARAMETERS

        USE density,    ONLY : derho,rho_cur,rho_new
        USE pressure,   ONLY : dep,p_cur,p_new
        USE u_velocity, ONLY : deu,u_cur,u_new
        USE v_velocity, ONLY : dev,v_cur,v_new
        USE energy,     ONLY : dee,e_cur,e_new

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ts
        
        REAL (KIND=8) :: rhox,ux,vx,px,ex
        REAL (KIND=8) :: rhoy,uy,vy,py,ey
        REAL (KIND=8) :: drho,du,dv,dp,de
        REAL (KIND=8) :: rho,u,v,p,e
        REAL (KIND=8) :: L1,L2,L3,L4

        INTEGER :: i,j

        rhox = 0.0d0
        ux   = 0.0d0
        vx   = 0.0d0
        px   = 0.0d0
        ex   = 0.0d0

        rhoy = 0.0d0
        uy   = 0.0d0
        vy   = 0.0d0
        py   = 0.0d0
        ey   = 0.0d0
        
        drho = 0.0d0
        du   = 0.0d0
        dv   = 0.0d0
        dp   = 0.0d0
        de   = 0.0d0

        rho = 0.0d0
        u   = 0.0d0
        v   = 0.0d0
        p   = 0.0d0
        e   = 0.0d0

        L1 = 0.0d0
        L2 = 0.0d0
        L3 = 0.0d0
        L4 = 0.0d0
        j = 1
        do i = 2,NX-1
        
        rhox = 0.5d0*(rho_new(i+1,1)-rho_new(i-1,1))/dx
        ux   = 0.5d0*(u_new(i+1,1)-u_new(i-1,1))/dx
        vx   = 0.5d0*(v_new(i+1,1)-v_new(i-1,1))/dx
        px   = 0.5d0*(p_new(i+1,1)-p_new(i-1,1))/dx
        ex   = 0.5d0*(e_new(i+1,1)-e_new(i-1,1))/dx

        rhoy = 0.5d0*(rho_new(i,2)-rho_new(i,1))/dy
        uy   = 0.5d0*(u_new(i,2)-u_new(i,1))/dy
        vy   = 0.5d0*(v_new(i,2)-v_new(i,1))/dy
        py   = 0.5d0*(p_new(i,2)-p_new(i,1))/dy
        ey   = 0.5d0*(e_new(i,2)-e_new(i,1))/dy

        c = sqrt(gama*p_new(i,j)/rho_new(i,j))

        if ((v_new(i,j)-c) >= 0.0d0) then
        L1 = 0.0d0
        else
        L1 = (v_new(i,j) - c)*(py - rho_new(i,j)*c*vy)
        end if

        if(v_new(i,j) >= 0.0d0) then
        L2 = 0.0d0
        L3 = 0.0d0
        else
        L2 = v_new(i,j)*(c*c*rhoy - py)
        L3 = v_new(i,j)*uy
        end if

        if((v_new(i,j)+c)>= 0.0d0) then
        L4 = 0.0d0
        else
        L4 = (v_new(i,j)+c)*(py+rho_new(i,j)*c*vy)
        end if



        drho = (L2+0.5d0*(L4+L1))/(c*c) + u_new(i,j)*rhox 
        drho = drho + rho_new(i,j)*ux

        dp = (L4 + L1)*0.5d0 + u_new(i,j)*px
        dp = dp + gama*p_new(i,j)*ux

        du = L3 + u_new(i,j)*ux + px/rho_new(i,j)

        dv = 0.5d0*(L4-L1)/(rho_new(i,j)*c)+u_new(i,j)*vx



        
        rho = rho_new(i,j) + drho * dt

        u = u_new(i,j) + du*dt

        v = v_new(i,j) + dv*dt

        p = p_new(i,j) + dp*dt

        e = p/(rho*(gama - 1.0D0))



        derho(ts,i) = 100.0d0*(rho_cur(i,j) - rho)/rho_cur(i,j)
        deu(ts,i)   = 100.0d0*(u_cur(i,j) - u)/u_cur(i,j)
        dev(ts,i)   = 100.0d0*(v_cur(i,j) - v)/v_cur(i,j)
        dep(ts,i)   = 100.0d0*(p_cur(i,j) - p)/p_cur(i,j)
        dee(ts,i)   = 100.0d0*(e_cur(i,j) - e)/e_cur(i,j)

        

        end do

        END SUBROUTINE
