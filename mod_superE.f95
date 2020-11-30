
      MODULE superE
!
! -----------------------------------------------------------------------------------------
!
      TYPE SuperElType
        REAL(8) a,b,c,e1,e2
        REAL(8) vos                  ! size of sphere enclosing the superellipsoid
        REAL(8) volume,mass,density
        REAL(8) moi(3,3)             ! moment of inertia tensor
        REAL(8) RM(3,3),RMT(3,3)     ! rotation matrix and its transpose
        REAL(8) ea(3) ! euler angles (phi,theta,psi)
        REAL(8) ep(4) ! euler parameters (e0,e1,e2,e3)
        REAL(8) r(3),rOld(3) ! location, and prev timestep location
        REAL(8) v(3),vOld(3)  ! velocity, and prev timestep velocity
        REAL(8) o(3)  ! angular velocity
        REAL(8) axis(3) ! (0,0,1) transformed from PRF into GRF
        INTEGER element ! mesh element in which particle is located
        INTEGER loctet  ! local tetrahedron
        INTEGER potentialElements(5) !in particle searc, potential found elements
        LOGICAL active ! particle is active or not
        INTEGER ierr ! particle error flag
        REAL(8) lambda ! b/a
        REAL(8) rdfx,rdfy,rdfz ! multiplication term in domega/dt eq.
        REAL(8) f,g,ksi,eta,hi !  elements of deformation rate tensor (f,g)
                               !  and the spin tensor (xi, eta, chi) in particle frame of reference
        REAL(8) AA,RR ! density ratios
        REAL(8) vs(3) ! settling velocity
        REAL(8) tau,St ! Particle response time, Stokes number
        REAL(8) taue,Ste ! Elliptic Particle response time, Elliptic Stokes number
        REAL(8) ResTprime(3) ! Resistance tensor prime (in PFR)
        INTEGER pid ! processor ID number
        INTEGER id ! my ID number

        !Mitja added
        INTEGER foundat !found ad which section (0 not found, 1 elem, 2 nei, 3 2nd-nei, 4 whole mesh)
        INTEGER nnf !number of not found occurences in a row
        INTEGER maxnnf !max number of not found, before deactivate particle
        !REAL(8) lastLoc(3) !location of particle when last found
        INTEGER atbound !if particle is found in boundary element

        !INTEGER cellHist(100000,2) !particle cell occupancy history        

      END TYPE

      CONTAINS

! -----------------------------------------------------------------------------------------
      SUBROUTINE seWriteBIN(lun,se)
!
!     Write single SE data to binary file
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se
      INTEGER lun

      WRITE(lun) se%id,se%element
      WRITE(lun) se%active,se%ierr
      WRITE(lun) se%a,se%b,se%c,se%e1,se%e2,se%density
      WRITE(lun) se%ep,se%r,se%v,se%o

      END SUBROUTINE

! -----------------------------------------------------------------------------------------
      SUBROUTINE seReadBIN(lun,se)
!
!     Write single SE data to binary file
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se
      INTEGER lun

      READ(lun) se%id,se%element
      READ(lun) se%active,se%ierr
      READ(lun) se%a,se%b,se%c,se%e1,se%e2,se%density
      READ(lun) se%ep,se%r,se%v,se%o

      END SUBROUTINE
! -----------------------------------------------------------------------------------------
      SUBROUTINE seCalVOS(se)
!
!     Calculate size of shpere enclosing the super ellipsoid
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se

!     an approximation (for low epsilon wrong)
      se%vos = MAX(se%a,se%b)
      se%vos = MAX(se%c,se%vos)

      se%vos=se%vos*SQRT(3.0D0) ! square is worst case scenarion (eps=0)


      END SUBROUTINE


! -----------------------------------------------------------------------------------------
      SUBROUTINE seCalRotationMatrix(se)
!
!     Forms the rotation matrix RM and its transpose RMtrans based on
!     Euler parameters in se%EP
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se

!     Rotation Matrix

      se%RM(1,1)=  se%ep(1)**2 + se%ep(2)**2 - se%ep(3)**2 - se%ep(4)**2
      se%RM(1,2)=  2.0D0*(se%ep(2)*se%ep(3) + se%ep(1)*se%ep(4))
      se%RM(1,3)=  2.0D0*(se%ep(2)*se%ep(4) - se%ep(1)*se%ep(3))

      se%RM(2,1)=  2.0D0*(se%ep(2)*se%ep(3) - se%ep(1)*se%ep(4))
      se%RM(2,2)=  se%ep(1)**2 - se%ep(2)**2 + se%ep(3)**2 - se%ep(4)**2
      se%RM(2,3)=  2.0D0*(se%ep(3)*se%ep(4) + se%ep(1)*se%ep(2))

      se%RM(3,1)=  2.0D0*(se%ep(2)*se%ep(4) + se%ep(1)*se%ep(3))
      se%RM(3,2)=  2.0D0*(se%ep(3)*se%ep(4) - se%ep(1)*se%ep(2))
      se%RM(3,3)=  se%ep(1)**2 - se%ep(2)**2 - se%ep(3)**2 + se%ep(4)**2

!     Transpose Rotation Matrix

      se%RMT(1,1)=  se%RM(1,1)
      se%RMT(1,2)=  se%RM(2,1)
      se%RMT(1,3)=  se%RM(3,1)

      se%RMT(2,1)=  se%RM(1,2)
      se%RMT(2,2)=  se%RM(2,2)
      se%RMT(2,3)=  se%RM(3,2)

      se%RMT(3,1)=  se%RM(1,3)
      se%RMT(3,2)=  se%RM(2,3)
      se%RMT(3,3)=  se%RM(3,3)

      END SUBROUTINE


! -----------------------------------------------------------------------------------------
      SUBROUTINE seOrientationAxis(se)
!
!     Based on Euler parameters (e0,e1,e2,e3) and rotation matrix
!     calculate particle orientation axis
!     Transforms (0,0,1) from PFR to GFR
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se
      REAL(8) cs(3)

!     In local (particle) coordinate system
      cs(1)=0.0D0
      cs(2)=0.0D0
      cs(3)=1.0D0

!     Transform to global cs
      se%axis=MATMUL(se%RMT,cs)

      END SUBROUTINE



! -----------------------------------------------------------------------------------------
      SUBROUTINE seMOI(se)
!
!     moment of inertia of superellipsoid (eq 2.67 in Jaklic : Segmentation and Recovery of Superquadrics)
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se

      se%MoI=0.0D0

      se%MoI(1,1) = 0.5D0 * se%a * se%b * se%c * se%e1 * se%e2 * ( &
     &              se%b*se%b*seGetBeta(1.5D0*se%e2,0.5D0*se%e2)* &
     &                        seGetBeta(0.5D0*se%e1,1.0D0+2.0D0*se%e1)+ &
     &        4.0D0*se%c*se%c*seGetBeta(0.5D0*se%e2,1.0D0+0.5D0*se%e2)* &
                             seGetBeta(1.5D0*se%e1,1.0D0+1.0D0*se%e1) )

      se%MoI(2,2) = 0.5D0 * se%a * se%b * se%c * se%e1 * se%e2 * ( &
     &              se%a*se%a*seGetBeta(1.5D0*se%e2,0.5D0*se%e2)* &
     &                        seGetBeta(0.5D0*se%e1,1.0D0+2.0D0*se%e1)+ &
     &        4.0D0*se%c*se%c*seGetBeta(0.5D0*se%e2,1.0D0+0.5D0*se%e2)* &
                             seGetBeta(1.5D0*se%e1,1.0D0+1.0D0*se%e1) )

      se%MoI(3,3) = 0.5D0 * se%a * se%b * se%c * se%e1 * se%e2 * &
     &              (se%a*se%a+se%b*se%b)* &
     &                        seGetBeta(1.5D0*se%e2,0.5D0*se%e2)* &
                             seGetBeta(0.5D0*se%e1,1.0D0+2.0D0*se%e1)


      END SUBROUTINE



! -----------------------------------------------------------------------------------------
      SUBROUTINE seSetVolume(se)
!
!     Calculate volume of superellipsoid (eq 2.59 in Jaklic : Segmentation and Recovery of Superquadrics)
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se

      se%volume=2.0D0*se%a*se%b*se%c*se%e1*se%e2* &
     &          seGetBeta(1.0D0+0.5D0*se%e1,se%e1)* &
               seGetBeta(0.5D0*se%e2,0.5D0*se%e2)

      END SUBROUTINE

! -----------------------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION seGetBeta(x,y)
!
!     Beta function is needed for volume and MoI calculations (eq 2.53 in Jaklic : Segmentation and Recovery of Superquadrics)
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) x,y

      seGetBeta = Gamma(x)*Gamma(y)/Gamma(x+y)

      END FUNCTION

! -----------------------------------------------------------------------------------------
      SUBROUTINE sePrint(se)
!
!     Print superellipsoid data on screen
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se

      print *,se%a,se%b,se%c,se%e1,se%e2
      print *,se%volume,se%mass,se%density
      print *,se%MoI(1,1)
      print *,se%MoI(2,2)
      print *,se%MoI(3,3)
      print *,se%ep

      END SUBROUTINE

! -----------------------------------------------------------------------------------------
      SUBROUTINE seTransfG2L(se,r)
!
!     Transform global point to local point
!     (IN)  r ! point in global (x,y,z) coordiante system
!     (OUT) r ! point in super ellipsoid local coordinate system (x',y',z')
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se

      REAL(8) r(3)

      r = MATMUL (se%RM,r-se%r)

      END SUBROUTINE
! -----------------------------------------------------------------------------------------
      SUBROUTINE seTransfL2G(se,r)
!
!     Transform global point to local point
!     (OUT)  r ! point in global (x,y,z) coordiante system
!     (IN) r ! point in super ellipsoid local coordinate system (x',y',z')
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se
      REAL(8) r(3)

      r = se%r + MATMUL(se%RMT,r)

      END SUBROUTINE


! -----------------------------------------------------------------------------------------
      SUBROUTINE seMapSphereToSuperE(se,t,r)
!
!     t = point on super ellipsoid (local ks)
!     r = point on sphere (|r|=1,c=(0,0,0)) (local ks)
!     (eq 2.10 in Jaklic : Segmentation and Recovery of Superquadrics)
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se
      REAL(8) t(3),r(3)
      REAL(8) eta,om,xy

      om=DATAN2(r(2),r(1))
      xy=SQRT(r(1)**2+r(2)**2)
      eta=DATAN2(r(3),xy)

      ! SIGN(A,B) returns the value of A with the sign of B
      t(1)=se%a * SIGN( (ABS(COS(eta)))**se%e1 , COS(eta) ) * SIGN( (ABS(COS(om)))**se%e2 , COS(om) )
      t(2)=se%b * SIGN( (ABS(COS(eta)))**se%e1 , COS(eta) ) * SIGN( (ABS(SIN(om)))**se%e2 , SIN(om) )
      t(3)=se%c * SIGN( (ABS(SIN(eta)))**se%e1 , SIN(eta) )

      END SUBROUTINE


! -----------------------------------------------------------------------------------------
      SUBROUTINE seEulerAngles2Parameters(se)
!
!     Transforms Euler angles to Euler parameters
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se

!     http://mathworld.wolfram.com/EulerParameters.html
!     3 = phi, 1= psi, 2= Theta
      se%ep(1) =   Cos(0.5D0*( se%ea(3) + se%ea(1)))*Cos(0.5D0*se%ea(2))
      se%ep(2) =   Cos(0.5D0*( se%ea(3) - se%ea(1)))*Sin(0.5D0*se%ea(2))
      se%ep(3) =   Sin(0.5D0*( se%ea(3) - se%ea(1)))*Sin(0.5D0*se%ea(2))
      se%ep(4) =   Sin(0.5D0*( se%ea(3) + se%ea(1)))*Cos(0.5D0*se%ea(2))

      END SUBROUTINE



! -----------------------------------------------------------------------------------------
      SUBROUTINE seEulerParameters2AnglesRotMat(se)
!
!     Transforms Euler parameters to Euler angles by using the rotation matrix
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      TYPE(SuperElType) se

!     gets euler angles via rotation matrix
      REAL(8) dvaPi
      dvaPi=8.0D0*DATAN(1.0D0)

!      zaradi napake pri eulerjevih parametrih lahko pride RM(3,3) malo nad 1 ali malo pod -1
      IF (se%RM(3,3).GT.1.01D0.OR.se%RM(3,3).LT.-1.01D0) PRINT *,"PAZI, EulerParameters2AnglesRotMat",se%RM(3,3)
      IF (se%RM(3,3).GT.1.0D0) se%RM(3,3)=1.0D0
      IF (se%RM(3,3).LT.-1.0D0) se%RM(3,3)=-1.0D0

!     Using the x-convention, the 3-1-3 extrinsic Euler angles
!     φ, θ and ψ (around the z-axis, x-axis and again the Z-axis) can be obtained as follows:
!     dela cudno za theta = ea(2) = 0.0

!     3 = phi, 1= psi, 2= Theta
      se%ea(1)=DATAN2(se%RM(1,3),se%RM(2,3))
      se%ea(2)=DACOS(se%RM(3,3))
      se%ea(3)=DATAN2(se%RM(3,1),-se%RM(3,2))

      IF (se%ea(1).LT.0.0D0) se%ea(1)=se%ea(1)+dvapi
      IF (se%ea(3).LT.0.0D0) se%ea(3)=se%ea(3)+dvapi

      END SUBROUTINE





! #############################################################################
      SUBROUTINE CalRotDynFactor(rdfx,rdfy,rdfz,a,lambda,f_visk,f_rho,f_L,p_rho,f_u0)
!
!     Calculates multiplication term in the domega/dt equations
!
! #############################################################################
      REAL(8) alfa,beta,gamal2,lambda,a,rdfx,rdfy,rdfz
      REAL(8) f_visk,f_rho,f_L,p_rho,f_u0


!     Calculates Gallily and Cohen parameters (note that gama is multiplied by lambda^2)
      CALL CalGalCoh(alfa,beta,gamal2,lambda)

!     x smer
      rdfx=20.0D0*f_visk*f_rho*f_L / ( a**2.0D0 * (beta+gamal2) * p_rho * f_u0 )    
!     y smer
      rdfy=20.0D0*f_visk*f_rho*f_L / ( a**2.0D0 * (alfa+gamal2) * p_rho * f_u0 )          
!     z smer
      rdfz=20.0D0*f_visk*f_rho*f_L / ( a**2.0D0 * (alfa+beta  ) * p_rho * f_u0 )          

      END SUBROUTINE


! #############################################################################
      SUBROUTINE CalGalCoh(alfa,beta,gama,lambda)
!
!     Calculates Gallily and Cohen parameters 
!     (note that gama is multiplied by lambda^2)
!
! #############################################################################
      REAL(8) alfa,beta,gama,lambda,l2

      IF (lambda.LT.1.0D0) THEN
        Print *,"Error in CalSt"
      ELSE IF (lambda.EQ.1.0D0) THEN
	  alfa=2.0D0/3.0D0
	  beta=2.0D0/3.0D0
	  gama=2.0D0/3.0D0
      ELSE
	  l2=lambda**2.0D0
	  alfa=l2/(l2-1.0D0)+lambda/(2.0D0*(l2-1.0D0)**1.5D0)*Log( (lambda-sqrt(l2-1.0D0))/(lambda+sqrt(l2-1.0D0)) )
	  beta=alfa
	  gama=-2.0D0*l2/(l2-1.0D0)-lambda*l2/((l2-1.0D0)**1.5D0)*Log( (lambda-sqrt(l2-1.0D0))/(lambda+sqrt(l2-1.0D0)) )
      END IF

      END SUBROUTINE


      END MODULE
