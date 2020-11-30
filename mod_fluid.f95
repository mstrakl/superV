!
! -----------------------------------------------------------------------------------------
!
      MODULE mFluid
!
! -----------------------------------------------------------------------------------------
!      
      TYPE fluidType
        INTEGER :: nelem,nnodes
       
        REAL(8), POINTER :: Un(:,:)    ! nnodes,3 :: velocity, nodal values
        REAL(8), POINTER :: Tn(:)    ! nnodes :: temperature, nodal values
        REAL(8), POINTER :: Vortn(:,:)    ! nnodes,3 :: vorticity, nodal values
        REAL(8), POINTER :: gradUxn(:,:)    ! nnodes,3 :: gradient of X velocity component, nodal values
        REAL(8), POINTER :: gradUyn(:,:)    ! nnodes,3 :: gradient of Y velocity component, nodal values
        REAL(8), POINTER :: gradUzn(:,:)    ! nnodes,3 :: gradient of Z velocity component, nodal values
        REAL(8), POINTER :: Ue(:,:)    ! nelem,3 :: velocity, element values
        REAL(8), POINTER :: Pn(:)    ! nnodes :: pressure, nodal values
        REAL(8), POINTER :: Pe(:)    ! nelem :: pressure, element values
        REAL(8), POINTER :: dvxdt(:),dvydt(:),dvzdt(:)

        INTEGER iUn,iUe,iPn,iPe,iVortn,iTn

      END TYPE fluidType 
!
! -----------------------------------------------------------------------------------------
!      
      TYPE fluidAtPartType
        REAL(8) vx,vy,vz  ! hitrost tekocine / fluid velocity
        REAL(8) wx,wy,wz  ! vrtincnost / vorticity
        REAL(8) gradVx(3),gradVy(3),gradVz(3) ! velocity gradients
        REAL(8) dvxdt,dvydt,dvzdt ! velocity time derivatives
      END TYPE
!
! -----------------------------------------------------------------------------------------
!

      END MODULE
