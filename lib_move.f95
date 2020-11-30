
! -----------------------------------------------------------------------------------------
      SUBROUTINE MoveParticle(part,fap,cpu)
!
!     Moves particles
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mFluid   
      USE logFile
      USE cpuTime      
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part

      REAL(8) c,cpu ! cpu time (this routine, total for this routine)
!
!     Init time measurement
!
      CALL cpStart(c)
!      
!     Update "Old" values      
!
      part%rOld=part%r 
      part%vOld=part%v
!
!     Move
!  
      IF (inp%MoveModel.EQ.0) THEN
        CALL MoveMasslessParticle(part,fap)
      ELSE IF (inp%MoveModel.EQ.1) THEN
        CALL MoveParticleRK4(part,fap)        
      ELSE IF (inp%MoveModel.EQ.2) THEN
        CALL MoveParticleEuler(part,fap)        
      ELSE
        CALL logIntWrite ("ERROR :: MoveParticle :: Unknown particle move model : ",inp%MoveModel)
        CALL StopProgram(1)
      END IF  
!
!     Stop time measurement
!
      CALL cpStop(c)
      cpu = cpu + c

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE MoveParticleRK4(part,fap)
!
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mFluid  
      USE counters 
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part
      INTEGER n

!     runge kutta variables
      REAL(8), ALLOCATABLE :: y(:),yout(:),dydx(:)

      n=13 ! stevilo enacb, number of equations to integrate by Runge-Kutta
      ALLOCATE (y(n),yout(n),dydx(n))  ! Runge - Kutta variables

      y(1)=part%r(1)  ! zacetni polozaj x, old position (x direction)
      y(2)=part%v(1) ! zacetna hitrost x, old velocity (x direction)
      y(3)=part%r(2) ! zacetni polozaj y old position (y direction)
      y(4)=part%v(2) ! zacetna hitrost y, old velocity (y direction)
      y(5)=part%r(3) ! zacetni polozaj z old position (z direction)
      y(6)=part%v(3) ! zacetna hitrost z, old velocity (z direction)
      y(7)=part%o(1) ! kotna hitrost v koordinatnem sistemu delca
      y(8)=part%o(2) ! kotna hitrost v koordinatnem sistemu delca
      y(9)=part%o(3) ! kotna hitrost v koordinatnem sistemu delca
      y(10)=part%ep(1) ! usmerjenost delca - Eulerjevi parametri
      y(11)=part%ep(2) ! usmerjenost delca - Eulerjevi parametri
      y(12)=part%ep(3) ! usmerjenost delca - Eulerjevi parametri
      y(13)=part%ep(4) ! usmerjenost delca - Eulerjevi parametri
      yout=y
!
!     Calculate derivatives of field functions
!
      CALL cal_dvdt(part,fap,y,dydx)
!	
!     Runge Kutta
!
      CALL rk4(y,dydx,n,cnt%rTime,inp%TimeStep,yout,part,fap)           
!
!     Store new particle data
!
      part%r(1)=yout(1)  
      part%v(1)=yout(2) 
      part%r(2)=yout(3)
      part%v(2)=yout(4)
      part%r(3)=yout(5)
      part%v(3)=yout(6)  
      part%o(1)=yout(7) ! kotna hitrost v koordinatnem sistemu delca
      part%o(2)=yout(8) ! kotna hitrost v koordinatnem sistemu delca
      part%o(3)=yout(9) ! kotna hitrost v koordinatnem sistemu delca
      part%ep(1)=yout(10) ! usmerjenost delca - Eulerjevi parametri
      part%ep(2)=yout(11) ! usmerjenost delca - Eulerjevi parametri
      part%ep(3)=yout(12) ! usmerjenost delca - Eulerjevi parametri
      part%ep(4)=yout(13) ! usmerjenost delca - Eulerjevi parametri

      DEALLOCATE (y,yout,dydx)

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE MoveParticleEuler(part,fap)
!
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mFluid  
      USE counters 
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part
      INTEGER n

!     runge kutta variables
      REAL(8), ALLOCATABLE :: y(:),yout(:),dydx(:)

      n=13 ! stevilo enacb, number of equations to integrate by Runge-Kutta
      ALLOCATE (y(n),yout(n),dydx(n))  ! Runge - Kutta variables

      y(1)=part%r(1)  ! zacetni polozaj x, old position (x direction)
      y(2)=part%v(1) ! zacetna hitrost x, old velocity (x direction)
      y(3)=part%r(2) ! zacetni polozaj y old position (y direction)
      y(4)=part%v(2) ! zacetna hitrost y, old velocity (y direction)
      y(5)=part%r(3) ! zacetni polozaj z old position (z direction)
      y(6)=part%v(3) ! zacetna hitrost z, old velocity (z direction)
      y(7)=part%o(1) ! kotna hitrost v koordinatnem sistemu delca
      y(8)=part%o(2) ! kotna hitrost v koordinatnem sistemu delca
      y(9)=part%o(3) ! kotna hitrost v koordinatnem sistemu delca
      y(10)=part%ep(1) ! usmerjenost delca - Eulerjevi parametri
      y(11)=part%ep(2) ! usmerjenost delca - Eulerjevi parametri
      y(12)=part%ep(3) ! usmerjenost delca - Eulerjevi parametri
      y(13)=part%ep(4) ! usmerjenost delca - Eulerjevi parametri
      yout=y
!
!     Calculate derivatives of field functions
!
      CALL cal_dvdt(part,fap,y,dydx)
!	
!     Runge Kutta
!          
      CALL EulerODE(y,dydx,n,inp%TimeStep,yout)
!
!     Store new particle data
!
      part%r(1)=yout(1)  
      part%v(1)=yout(2) 
      part%r(2)=yout(3)
      part%v(2)=yout(4)
      part%r(3)=yout(5)
      part%v(3)=yout(6)  
      part%o(1)=yout(7) ! kotna hitrost v koordinatnem sistemu delca
      part%o(2)=yout(8) ! kotna hitrost v koordinatnem sistemu delca
      part%o(3)=yout(9) ! kotna hitrost v koordinatnem sistemu delca
      part%ep(1)=yout(10) ! usmerjenost delca - Eulerjevi parametri
      part%ep(2)=yout(11) ! usmerjenost delca - Eulerjevi parametri
      part%ep(3)=yout(12) ! usmerjenost delca - Eulerjevi parametri
      part%ep(4)=yout(13) ! usmerjenost delca - Eulerjevi parametri

      DEALLOCATE (y,yout,dydx)

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE Acc_Gravity(part,dvdt)
!
!     Gravity and bouyancy acceleration of the particle
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mFluid   
      IMPLICIT NONE

      INTEGER i
      TYPE(SuperElType) part
      REAL(8) dvdt(3)

      DO i=1,3     
        dvdt(i)=part%AA/part%st*part%vs(i)
      END DO

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE Acc_StokesDrag(part,fap,cpv,dvdt)
!
!     Stokes drag
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mFluid   
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part
      REAL(8) dvdt(3)
      REAL(8) cpv(3)
   
      dvdt(1)=part%AA/part%st*(fap%vx-cpv(1))
      dvdt(2)=part%AA/part%st*(fap%vy-cpv(2))
      dvdt(3)=part%AA/part%st*(fap%vz-cpv(3))

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE Acc_EllipticDrag(part,fap,cpv,dvdt)
!
!     Elliptic drag (Stokes part odštet)
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mFluid   
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part
      REAL(8) dvdt(3)
      REAL(8) cpv(3),ResT(3,3)


!     calculate resistance in inertial frame of reference divided by (6*lambda)
      CALL CalResistanceTensor(ResT,part)

!     take care of separation of spherical and elliptical part (so we can use the code as spherical)
      ResT(1,1)=ResT(1,1)-1.0D0
      ResT(2,2)=ResT(2,2)-1.0D0
      ResT(3,3)=ResT(3,3)-1.0D0

      dvdt(1) = part%AA/part%st*(ResT(1,1)*(fap%vx-cpv(1))+ResT(1,2)*(fap%vy-cpv(2))+ResT(1,3)*(fap%vz-cpv(3)))
      dvdt(2) = part%AA/part%st*(ResT(2,1)*(fap%vx-cpv(1))+ResT(2,2)*(fap%vy-cpv(2))+ResT(2,3)*(fap%vz-cpv(3)))
      dvdt(3) = part%AA/part%st*(ResT(3,1)*(fap%vx-cpv(1))+ResT(3,2)*(fap%vy-cpv(2))+ResT(3,3)*(fap%vz-cpv(3)))

      END




! -----------------------------------------------------------------------------------------
      SUBROUTINE CalResistanceTensor(ResT,part)
!
!     Calculates resistance tensor REST  divided by (6*lambda) in inertial frame of reference
!     ResTprime - resistance tensor in particle frame of reference
!
! -----------------------------------------------------------------------------------------
      USE superE
      IMPLICIT NONE

      TYPE(SuperElType) part

      REAL(8) ResT(3,3)
      REAL(8), ALLOCATABLE :: tmp(:,:)
   	  
      INTEGER i,j

!     Calculate rotation matrix
      CALL seCalRotationMatrix(part)  ! v principu to izracuna ze getFFap
!     transform to inertial coor. system
      ALLOCATE (tmp(3,3))
      DO j=1,3
	    DO i=1,3
            tmp(i,j) = part%RM(i,j)*part%ResTprime(i)/6.0D0/part%lambda
  	    END DO
      END DO

      ResT=MATMUL(part%RMT,tmp)

      DEALLOCATE (tmp)

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE Acc_AmPc(part,fap,cpv,dvdt)
!
!     Added mass & pressure correction
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mFluid   
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part
      REAL(8) dvdt(3),cpv(3)


      dvdt(1)=1.5D0*part%RR*fap%dvxdt+ &
     &        part%RR*( &
     &               (fap%vx+0.5D0*cpv(1))*fap%gradVx(1)+ &
     &               (fap%vy+0.5D0*cpv(2))*fap%gradVx(2)+ &
                    (fap%vz+0.5D0*cpv(3))*fap%gradVx(3) )
      dvdt(2)=1.5D0*part%RR*fap%dvydt+ &
     &        part%RR*( &
     &               (fap%vx+0.5D0*cpv(1))*fap%gradVy(1)+ &
     &               (fap%vy+0.5D0*cpv(2))*fap%gradVy(2)+ &
                    (fap%vz+0.5D0*cpv(3))*fap%gradVy(3) )
      dvdt(3)=1.5D0*part%RR*fap%dvzdt+ &
     &        part%RR*( &
     &               (fap%vx+0.5D0*cpv(1))*fap%gradVz(1)+ &
     &               (fap%vy+0.5D0*cpv(2))*fap%gradVz(2)+ &
                    (fap%vz+0.5D0*cpv(3))*fap%gradVz(3) )

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE cal_dvdt(part,fap,y,dydx)
!
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mFluid   
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part

      REAL(8) dvdt(3),dydx(13),dodt(3),dedt(4),dodtmax
      REAL(8) cpv(3),y(13)  ! current particle velocity, to je tista, ki jo 
                            ! runge kutta popravlja, zacetna je v fap
                            ! ko je runge kutta skozi, potem fap%vx=cpv(1)
      REAL(8) Gravity(3),DragS(3),DragE(3),AmPC(3) ! accelerations
      REAL(8) l1,l2 ! values based on lambda
      REAL(8) ox,oy,oz,e0,e1,e2,e3
      INTEGER i


!     da ne pomesam
      cpv(1)=y(2) ! vx
      cpv(2)=y(4) ! vy
      cpv(3)=y(6) ! vz
      ox=y(7) ! ox
      oy=y(8) ! oy
      oz=y(9) ! oz
      e0=y(10) ! Euler parameter
      e1=y(11) ! Euler parameter
      e2=y(12) ! Euler parameter
      e3=y(13) ! Euler parameter

!
!     Accelerations
!
      DO i=1,3     
        Gravity(i)=0.0D0
        DragS(i)=0.0D0
        DragE(i)=0.0D0
        AmPC(i)=0.0D0
      END DO
      IF (inp%fm_Gravity.GT.0) CALL Acc_Gravity(part,Gravity)
      IF (inp%fm_StokesDrag.GT.0) CALL Acc_StokesDrag(part,fap,cpv,DragS)
      IF (inp%fm_EllipticDrag.GT.0) CALL Acc_EllipticDrag(part,fap,cpv,DragE)
      IF (inp%fm_AmPc.GT.0) CALL Acc_AmPc(part,fap,cpv,AmPC)
      DO i=1,3     
        dvdt(i)=Gravity(i)+DragS(i)+DragE(i)+AmPC(i)
      END DO

!
!     Particle angular velocity
!
      l1=(part%lambda**2.0D0-1.0D0) / ( 1.0D0 + part%lambda**2 )
!     smer x
      dodt(1)= oy*oz*l1+part%rdfx * (-l1*part%f + part%ksi - ox )
!     smer y
      dodt(2)=-oz*ox*l1+part%rdfy * ( l1*part%g + part%eta - oy )
!     smer z
      dodt(3)=          part%rdfz * (                 part%hi  - oz )

!
!     Sanity check :: angular velocity
!
      dodtmax = 1.0E+06
      DO i=1,3
            IF (dodt(i).GT.dodtmax) THEN
                  dodt(i) = dodtmax
                  !PRINT *, "Error:: dodt(",i,") > dodtmax! Correcting!"
            END IF
      END DO
!
!     Particle Rotation - Euler parameters
!
      dedt(1)=0.5D0 * ( -e1*ox - e2*oy - e3*oz )
      dedt(2)=0.5D0 * ( +e0*ox - e3*oy + e2*oz )
      dedt(3)=0.5D0 * ( +e3*ox + e0*oy - e1*oz )
      dedt(4)=0.5D0 * ( -e2*ox + e1*oy + e0*oz )
!
!     Copy all r.h.s. to Runge-Kutta field
!
      dydx(1)=y(2)    ! dx / dt = vx
      dydx(2)=dvdt(1) ! dvx / dt = ax
      dydx(3)=y(4)    ! dy / dt = vy
      dydx(4)=dvdt(2) ! dvy / dt = ay
      dydx(5)=y(6)    ! dz / dt = vz
      dydx(6)=dvdt(3) ! dvz / dt = az
      dydx(7)=dodt(1) ! dox/dt        
      dydx(8)=dodt(2) ! doy/dt
      dydx(9)=dodt(3) ! doz/dt
      dydx(10)=dedt(1) ! de0/dt
      dydx(11)=dedt(2) ! de0/dt
      dydx(12)=dedt(3) ! de0/dt
      dydx(13)=dedt(4) ! de0/dt

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE MoveMasslessParticle(part,fap)
!
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mFluid   
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part

      part%r(1)=part%r(1)+fap%vx*inp%TimeStep
      part%r(2)=part%r(2)+fap%vy*inp%TimeStep
      part%r(3)=part%r(3)+fap%vz*inp%TimeStep 

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE EulerODE(y,dydx,n,h,yout)
!
!     Euler ODE solver, 2th order accurate
!
! -----------------------------------------------------------------------------------------      
      USE superE
      USE mFluid   
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part

      INTEGER n,i
      REAL(8) h,dydx(n),y(n),yout(n)

      DO i=1,n
        yout(i)=y(i)+h*dydx(i)
      END DO

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE rk4(y,dydx,n,x,h,yout,part,fap)
!
!     Runge Kutta ODE solver, 4th order accurate
!
! -----------------------------------------------------------------------------------------      
      USE superE
      USE mFluid   
      IMPLICIT NONE

      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) part

      INTEGER n,i
      REAL*8 h,x,dydx(n),y(n),yout(n)
      REAL*8 h6,hh,xh,dym(n),dyt(n),yt(n)

      hh=h*0.5
      h6=h/6.
      xh=x+hh
      do 11 i=1,n
        yt(i)=y(i)+hh*dydx(i)
11    continue
      CALL cal_dvdt(part,fap,yt,dyt)
      do 12 i=1,n
        yt(i)=y(i)+hh*dyt(i)
12    continue
      CALL cal_dvdt(part,fap,yt,dym)
      do 13 i=1,n
        yt(i)=y(i)+h*dym(i)
        dym(i)=dyt(i)+dym(i)
13    continue
      CALL cal_dvdt(part,fap,yt,dyt)
      do 14 i=1,n
        yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.0D0*dym(i))
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.    


