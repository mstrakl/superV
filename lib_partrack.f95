! -----------------------------------------------------------------------------------------
      SUBROUTINE GetFFFap(m,p,fluid,fap,iT,ierr)
!
!     Interpolate flow field data to particle position
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      USE mFluid
      USE cpuTime
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) :: p ! particle
      TYPE(meshType) :: m  ! mesh data structure    
      TYPE(fluidType) :: fluid ! fluid flow fields
      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position

      INTEGER i,ie,ifound,ierr,iT
      REAL(8) cc,c ! cpu time (this routine, total for this routine)
      REAL(8) pr(3),xez(3)

      REAL(8) debugout(100,3)

      debugout = 0.0D0

      ierr = 0
!
!     Init time measurement
!
      CALL cpStart(c)
!
!     Find particle mesh element (only actual timeloop, when iT > 0)
!
      IF (iT.GT.1) THEN

        CALL cpStart(cc)
        !CALL FindParticleInUnstructuredMesh2(m,p,iT,ierr) 
        CALL FindParticleInUnstructuredMesh4(m,p,iT,ierr,debugout) 
        CALL cpStop(cc)
        cput%meas(4)=cput%meas(4) + cc

        IF (ierr.EQ.1) THEN 
        
          PRINT *, "!--------------------------------------------!"
          PRINT *, "GetFFFap:: After Find Particle ierr=",ierr
          WRITE(*,'(A10,I6)') "iT=",iT
          WRITE(*,'(A10,I6)') "Part id=",p%id
          WRITE(*,'(A10,3F10.4)') "Part r=",p%r
          WRITE(*,'(A10,I6)') "Element=",p%element
          PRINT *, "!--------------------------------------------!"
          !CALL BreakPoint()   
            
        END IF

        IF (ierr.NE.0) GOTO 10
      
      END IF
!
!     Interpolate within element
!
      CALL cpStart(cc)   

      pr = p%r
      ie = p%element
      
      CALL FindWithinElementTetDec(m,pr,ie,ifound,ierr)
      
      !Update if changed within FindWithinElementTetDec
      p%element = ie

      IF (ierr.EQ.0) THEN !found

        p%nnf = 0

      ELSE 

        IF (iT.GT.1 .AND. ierr.EQ.1) THEN

          PRINT *, "!--------------------------------------------!"
          PRINT *, "GetFFFap:: After In Element Find ierr=",ierr
          WRITE(*,'(A10,I6)') "iT=",iT
          WRITE(*,'(A10,I6)') "Part id=",p%id
          WRITE(*,'(A10,3F10.4)') "Part r=",p%r
          WRITE(*,'(A10,I6)') "Element=",ie
          PRINT *, "!--------------------------------------------!"
  
          !CALL WriteVTKTriangulatedElement('trielem',iT,m,ie)
          CALL WriteVTKPoint('diag/pntCC',iT,m%e(ie)%xc)
          CALL WriteVTKPoint('diag/pntA',iT,debugout(1,:))
          CALL WriteVTKPoint('diag/pntB',iT,debugout(2,:))

          DO i = 10,100

            IF (INT(debugout(i,1)).GT.0) THEN
              ie = INT(debugout(i,1))
              CALL WriteVTKTriangulatedElement('diag/trielem',ie,m,ie)
            END IF

          END DO
  
          p%nnf = p%nnf + 1
          !CALL BreakPoint()

        END IF
        
        GOTO 10     
          
      END IF

      CALL InterpolateFFInElementNew(m,p,pr,ie,ifound,fluid,fap)

      CALL cpStop(cc)
      cput%meas(6)=cput%meas(6) + cc  

      !CALL WriteParticleCellHistory(p,iT)
       
!
!     Calculate elements of strain and rotation tensor in particle frame of reference
!   
      CALL CalVGTpfr(fap,p)
!
!     Stop time measurement
!
10    CONTINUE
      CALL cpStop(c)
      cput%meas(3)=cput%meas(3) + c

      END

!! -----------------------------------------------------------------------------------------
!      SUBROUTINE WriteParticleCellHistory(p,iT)
!!
!!      
!! -----------------------------------------------------------------------------------------
!
!      USE superE
!      IMPLICIT NONE
!
!      TYPE(SuperElType) :: p ! particle
!
!      INTEGER :: iT
!      INTEGER :: i,ilast
!
!
!      IF (iT.LE.1) THEN
!        p%cellHist = -1
!
!        p%cellHist(1,1) = p%element
!        p%cellHist(1,2) = iT
!
!        RETURN
!      END IF
!
!      i=1
!      DO,WHILE ( p%cellHist(i,1).NE.-1 )
!        i = i + 1
!      END DO
!      ilast = i
!
!      IF (p%element.NE.p%cellHist(ilast-1,1)) THEN
!        
!        p%cellHist(ilast,1) = p%element
!        p%cellHist(ilast,2) = iT
!
!      END IF
!      
!      
!      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE interpolationTesting(m,p,fluid,fap,iT,ierr)
!
!     Interpolate flow field data to particle position
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      USE mFluid
      USE cpuTime
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) :: p ! particle
      TYPE(meshType) :: m  ! mesh data structure    
      TYPE(fluidType) :: fluid ! fluid flow fields
      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position

      INTEGER ierr,iT

      CHARACTER(255) fname

      INTEGER i,j,k,l,n,ndiv,iom

      REAL(8) xmin,xmax,ymin,ymax,zmin,zmax
      REAL(8), ALLOCATABLE :: vv(:,:), outMatrix(:,:)

      REAL(8) minabsolute,varAna(3),varInt(3)
      REAL(8) v(3),uint(3),uana(3),err1,err2,err3 !,xval,yval,zval,uval
      REAL(8) gradUzAna(3), gradUzInt(3)

      REAL(8) fac,pi

      REAL(8) cput1,cput2

      REAL(8) xez(3)
     
      p%element = 21

      p%r(1) = 0.0065D0
      p%r(2) = 0.01D0
      !p%r(3) = 0.089D0
      p%r(3) = 0.093D0

      !CALL GetFFFap(m,p,fluid,fap,iT,ierr)

      IF (m%e(p%element)%type.NE.42) THEN

        CALL GetXiEtaZeta(m,m%e(p%element),p%r,xez,ierr)

        IF (ierr.NE.0) THEN
          PRINT *,"Lost particle at XiEtaZeta" 
        END IF 
        WRITE(*,'(A4,3G12.4)') "XEZ=",xez

      ELSE

        xez=p%r

      END IF  
!
!     Interpolate within element
!

      !CALL InterpolateFFInElementNew(m,m%e(p%element),p%element,xez,fluid,fap,ierr)

      uint(1) = fap%vx
      uint(2) = fap%vy
      uint(3) = fap%vz    

      WRITE(*,'(A5,3F16.6)') "uint=",uint

      PRINT *, "Stopping in int test"
      CALL stopProgram(1)


      END      



! -----------------------------------------------------------------------------------------
      SUBROUTINE CalVGTpfr(fap,p)
!
!     Calculate elements of deformation rate tensor (f,g)
!     and the spin tensor (xi, eta, chi) in particle frame of reference
!     gradVx,gradVy,gradVz - velocity gradient in inertial frame of reference
!     vgt - velocity gradient tensor in particle frame of reference
!      
! -----------------------------------------------------------------------------------------
      USE SuperE
      USE mFluid
      IMPLICIT NONE

      TYPE(SuperElType) p
      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position

      REAL(8), ALLOCATABLE :: tmp(:,:),vgt(:,:),gradV(:,:)

      ALLOCATE (tmp(3,3),vgt(3,3))
      ALLOCATE (gradV(3,3))
!
!     Calculate rotation matrix and its transpose
!
      CALL seCalRotationMatrix(p)

!     Calculate velocity gradient tensor in particle frame of reference
      gradV(1,1)=fap%gradVx(1)
      gradV(1,2)=fap%gradVx(2)
      gradV(1,3)=fap%gradVx(3)

      gradV(2,1)=fap%gradVy(1)
      gradV(2,2)=fap%gradVy(2)
      gradV(2,3)=fap%gradVy(3)

      gradV(3,1)=fap%gradVz(1)
      gradV(3,2)=fap%gradVz(2)
      gradV(3,3)=fap%gradVz(3)

!     transform from inertial to particle frame of reference
      tmp=MATMUL(gradV,p%RMT)
      vgt=MATMUL(p%RM,tmp)


!     Get strain and rotation tensor elements from velocity gradient tensor in particle frame of reference
      CALL CalTensorElements(vgt,p%f,p%g,p%ksi,p%eta,p%hi)


      DEALLOCATE(gradV,tmp,vgt) 
	
      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE CalTensorElements(vgt,f,g,ksi,eta,hi)
!
!     Form elements of deformation rate tensor (f,g)
!     and elements of the spin tensor (xi, eta, chi)
!     vgt - velocity gradient tensor in particle frame of reference
!
! -----------------------------------------------------------------------------------------
      REAL(8) vgt(3,3),f,g,ksi,eta,hi

      f  =0.5D0*( vgt(3,2) + vgt(2,3) )
      g  =0.5D0*( vgt(1,3) + vgt(3,1) )
      ksi=0.5D0*( vgt(3,2) - vgt(2,3) )
      eta=0.5D0*( vgt(1,3) - vgt(3,1) )
      hi =0.5D0*( vgt(1,2) - vgt(2,1) )

      END

! -----------------------------------------------------------------------------------------
      
      SUBROUTINE FindWithinElementTetDec(m,r,ie,ifound,ierr)
!
!     Search for particle within element, by decomposing in to tetrahedrons
!        
! -----------------------------------------------------------------------------------------
      
      USE mesh
      USE superE
      USE mFluid
      USE logFile
      USE counters

      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure

      INTEGER ie, ifound, ierr
      INTEGER i,j
      INTEGER iface,newie,newie2
      REAL(8) r(3)

      ierr = 1

!-----------------------------------------------------------      
!     Find local tet, which occupies the particle
!-----------------------------------------------------------

      ifound = findInTet(r,ie)

      IF (ifound.GT.0) GOTO 10

!-----------------------------------------------------------      
!     If not found, check neighbour elements
!-----------------------------------------------------------

      DO i = 1,size(m%e(ie)%faceCon)

        iface = m%e(ie)%faceCon(i)

        IF (iface.GT.0) THEN
          newie = m%f(abs(iface))%neighbour
        ELSE
          newie = m%f(abs(iface))%owner
        END IF

        ifound = findInTet(r,newie)

        IF (ifound.GT.0) THEN
          ie = newie
          GOTO 10
        END IF

      END DO

!-----------------------------------------------------------      
!     If not found, check second level neighbours
!-----------------------------------------------------------

      DO i = 1,size(m%e(ie)%faceCon)

        iface = m%e(ie)%faceCon(i)
        IF (iface.GT.0) THEN
          newie = m%f(abs(iface))%neighbour
        ELSE
          newie = m%f(abs(iface))%owner
        END IF

        DO j = 1,size(m%e(newie)%faceCon)

          iface = m%e(ie)%faceCon(i)
          IF (iface.GT.0) THEN
            newie2 = m%f(abs(iface))%neighbour
          ELSE
            newie2 = m%f(abs(iface))%owner
          END IF

          ifound = findInTet(r,newie2)

          IF (ifound.GT.0) THEN
            ie = newie2
            GOTO 10
          END IF

        END DO

      END DO

      !WRITE(*,'(A30,I6)') "For element :: ie", ie
      !WRITE(*,'(A4,3F12.6)') "r=", r
      !WRITE(*,'(A4,3F12.6)') "ec=", e%xc
      !PRINT *, "Local tet find results:"
      !WRITE(*,'(A8,I6)') "ifound=", ifound
      !WRITE(*,'(A8,I6)') "nfound=", k

      RETURN


10    CONTINUE
     
      ierr = 0

      !check if not boundary
      IF (ie.LT.1) THEN
        ierr = 2
        PRINT *, "Boundary in interpolate, thats strange"
        CALL BreakPoint()
      END IF


      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      INTEGER FUNCTION findInTet(r,ie)

        IMPLICIT NONE
        INTEGER, INTENT(IN)   :: ie
        REAL(8), INTENT(IN)   :: r(3)

        INTEGER i,ntet,ierr
        REAL(8) res,fCenters(4,3),fNormals(4,3)

        findInTet = 0

        IF (ie.GT.0) THEN

          ntet = size (m%e(ie)%tet)

          DO i = 1,ntet
    
            fCenters = m%e(ie)%tet(i)%centers
            fNormals = m%e(ie)%tet(i)%normals
    
            !PRINT *, "Searching loop:", i,"--------------------"
            !WRITE(*,'(A10,3F14.6)') "center 1=", fCenters(1,:)
            !WRITE(*,'(A10,3F14.6)') "center 2=", fCenters(2,:)
            !WRITE(*,'(A10,3F14.6)') "center 3=", fCenters(3,:)
            !WRITE(*,'(A10,3F14.6)') "center 4=", fCenters(4,:)
            !WRITE(*,'(A10,3F14.6)') "normal 1=", fNormals(1,:)
            !WRITE(*,'(A10,3F14.6)') "normal 2=", fNormals(2,:)
            !WRITE(*,'(A10,3F14.6)') "normal 3=", fNormals(3,:)
            !WRITE(*,'(A10,3F14.6)') "normal 4=", fNormals(4,:)
    
            CALL PartInElement2(r,fCenters,fNormals,4,res,ierr)
    
            IF (ierr.EQ.0) THEN
              findInTet = i 
              RETURN
            END IF
    
          END DO

        END IF

      END FUNCTION findInTet


      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE InterpolateFFInElementNew(m,p,r,ie,ifound,fluid,fap)
!
!     Interpolate flow field to particle position
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      USE mFluid
      USE logFile
      USE counters

      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure
      TYPE(fluidType) :: fluid ! fluid flow fields
      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) :: p ! particle

      INTEGER ie, dummierr
      REAL(8) r(3),xez(3)

      INTEGER i,j,ifound,inode
      INTEGER loctcon(3)

      REAL(8) x1(3),x2(3),x3(3),x4(3)
      REAL(8) fi(4)

      REAL(8),ALLOCATABLE :: uin1(:),uin2(:),uin3(:)

!-----------------------------------------------------------      
!     Get local tet con, xez and shape functions
!----------------------------------------------------------- 

      DO i = 1,m%e(ie)%nVertex

        IF  (m%e(ie)%nodecon(i).EQ.m%e(ie)%tet(ifound)%con(1) ) THEN
          loctcon(1) = i
        ELSE IF  (m%e(ie)%nodecon(i).EQ.m%e(ie)%tet(ifound)%con(2) ) THEN
          loctcon(2) = i
        ELSE IF  (m%e(ie)%nodecon(i).EQ.m%e(ie)%tet(ifound)%con(3) ) THEN
          loctcon(3) = i
        END IF 

      END DO

      x1 = m%x( m%e(ie)%tet(ifound)%con(1),: )
      x2 = m%x( m%e(ie)%tet(ifound)%con(2),: )
      x3 = m%x( m%e(ie)%tet(ifound)%con(3),: )
      x4 = m%e(ie)%xc

      CALL tet_kks2lks(r,x1,x2,x3,x4,xez,dummierr)

      CALL tet_shapef(xez,fi)

!-----------------------------------------------------------      
!     Interpolate in chosen tet with shape functions
!----------------------------------------------------------- 

      !PRINT *, "nvert=",m%e(ie)%nVertex
      !PRINT *, "con=",m%e(ie)%tet(ifound)%con
      !PRINT *, "lcon=",loctcon
      !PRINT *, "fi=",fi

      ALLOCATE( uin1(m%e(ie)%nVertex) )
      ALLOCATE( uin2(m%e(ie)%nVertex) )
      ALLOCATE( uin3(m%e(ie)%nVertex) )

      !
      !interpolate :: v
      !
      DO i = 1,m%e(ie)%nVertex
        inode = m%e(ie)%nodeCon(i)
        uin1(i) = fluid%Un(inode,1)
        uin2(i) = fluid%Un(inode,2)
        uin3(i) = fluid%Un(inode,3)
      END DO

      fap%vx = interpolate(m%e(ie)%nVertex, loctcon, uin1, fi)
      fap%vy = interpolate(m%e(ie)%nVertex, loctcon, uin2, fi)
      fap%vz = interpolate(m%e(ie)%nVertex, loctcon, uin3, fi)

      !
      !interpolate :: w
      !
      DO i = 1,m%e(ie)%nVertex
        inode = m%e(ie)%nodeCon(i)
        uin1(i) = fluid%Vortn(inode,1)
        uin2(i) = fluid%Vortn(inode,2)
        uin3(i) = fluid%Vortn(inode,3)
      END DO

      fap%wx = interpolate(m%e(ie)%nVertex, loctcon, uin1, fi)
      fap%wy = interpolate(m%e(ie)%nVertex, loctcon, uin2, fi)
      fap%wz = interpolate(m%e(ie)%nVertex, loctcon, uin3, fi)

      !
      !interpolate :: gradV
      !
      DO j = 1,3

        DO i = 1,m%e(ie)%nVertex
          inode = m%e(ie)%nodeCon(i)
          uin1(i) = fluid%gradUxn(inode,j)
          uin2(i) = fluid%gradUyn(inode,j)
          uin3(i) = fluid%gradUzn(inode,j)
        END DO

        fap%gradVx(j) = interpolate(m%e(ie)%nVertex, loctcon, uin1, fi)
        fap%gradVy(j) = interpolate(m%e(ie)%nVertex, loctcon, uin2, fi)
        fap%gradVz(j) = interpolate(m%e(ie)%nVertex, loctcon, uin3, fi)

      END DO

      !
      !interpolate :: dvdt
      !
      DO i = 1,m%e(ie)%nVertex
        inode = m%e(ie)%nodeCon(i)
        uin1(i) = fluid%dvxdt(inode)
        uin2(i) = fluid%dvydt(inode)
        uin3(i) = fluid%dvzdt(inode)
      END DO

      fap%dvxdt = interpolate(m%e(ie)%nVertex, loctcon, uin1, fi)
      fap%dvydt = interpolate(m%e(ie)%nVertex, loctcon, uin2, fi)
      fap%dvzdt = interpolate(m%e(ie)%nVertex, loctcon, uin3, fi)


      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      FUNCTION interpolate(np,con,uin,fi) RESULT(uout)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: np
        INTEGER, INTENT(IN) :: con(3)
        REAL(8), INTENT(IN) :: uin(np)
        REAL(8), INTENT(IN) :: fi(4)
        
        INTEGER i,inode
        REAL(8) ucent
        REAL(8) uout

        !get element center value
        ucent = 0.0D0
        DO i = 1,np
          ucent = ucent + uin(i)
        END DO

        ucent = ucent / np

        uout = 0.0D0
        DO i = 1,3
          inode = con(i)
          uout = uout + fi(i) * uin(inode)
        END DO

        !i=4 - center value
        uout = uout + fi(4) * ucent


    
      END FUNCTION interpolate

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE InterpolateFFInElement(m,p,e,ie,r,fluid,fap,ierr)
!
!     Interpolate flow field to particle position
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      USE mFluid
      USE logFile
      USE counters

      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure
      TYPE(elementType) :: e  ! mesh element
      TYPE(fluidType) :: fluid ! fluid flow fields
      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position
      TYPE(SuperElType) :: p ! particle

      INTEGER ie, ierr, dummierr
      REAL(8) r(3),xez(3),res,minres

      INTEGER i,j,k,ntet,ifound,inode,iminres
      INTEGER loctcon(3)

      REAL(8) fCenters(4,3), fNormals(4,3)
      REAL(8) x1(3),x2(3),x3(3),x4(3)
      REAL(8) fi(4)

      REAL(8) uin1(e%nVertex)
      REAL(8) uin2(e%nVertex)
      REAL(8) uin3(e%nVertex)

      ierr = 0

!-----------------------------------------------------------      
!     Find local tet, which occupies the particle
!-----------------------------------------------------------

      !init
      minres = 1E+10

      ntet = size (e%tet)

      !WRITE(*,'(A30,I6)') "Interpolating :: ie", ie
      !WRITE(*,'(A4,3F12.6)') "r=", r
      !WRITE(*,'(A4,3F12.6)') "ec=", e%xc

      k = 0
      ifound = 0
      DO i = 1,ntet

        fCenters = e%tet(i)%centers
        fNormals = e%tet(i)%normals

        !PRINT *, "Searching loop:", i,"--------------------"
        !WRITE(*,'(A10,3F14.6)') "center 1=", fCenters(1,:)
        !WRITE(*,'(A10,3F14.6)') "center 2=", fCenters(2,:)
        !WRITE(*,'(A10,3F14.6)') "center 3=", fCenters(3,:)
        !WRITE(*,'(A10,3F14.6)') "center 4=", fCenters(4,:)
        !WRITE(*,'(A10,3F14.6)') "normal 1=", fNormals(1,:)
        !WRITE(*,'(A10,3F14.6)') "normal 2=", fNormals(2,:)
        !WRITE(*,'(A10,3F14.6)') "normal 3=", fNormals(3,:)
        !WRITE(*,'(A10,3F14.6)') "normal 4=", fNormals(4,:)

        CALL PartInElement2(r,fCenters,fNormals,4,res,ierr)

        !save minres, for missed elements
        IF (res.LT.minres) THEN
          minres = res
          iminres = i
        END IF

        IF (ierr.EQ.0) THEN
          ifound = i
          k = k + 1
        END IF

      END DO

      !important leave this
      ierr = 0

      !PRINT *, "Local tet find results:"
      !WRITE(*,'(A8,I6)') "ifound=", ifound
      !WRITE(*,'(A8,I6)') "nfound=", k

!-----------------------------------------------------------      
!      If not found, deal with this shit somehow
!-----------------------------------------------------------

      IF (k.EQ.0) THEN !not found

        ifound = iminres
        cnt%inthowfound(2) = cnt%inthowfound(2)+1
        ierr = 1

      ELSE IF (k.EQ.1) THEN !found

        cnt%inthowfound(1) = cnt%inthowfound(1)+1

      ELSE IF (k.GT.1) THEN !multi found

        !happens when particles are very close
        !to the tet-element surface. It's not really
        !an error, since the particle is so close
        !to both element, that it doesn't matter which
        !element is used for the interpolation

        cnt%inthowfound(3) = cnt%inthowfound(3)+1

        !CALL WriteVTKTriangulatedElement('trielement.vtk',m,p%element)
        !CALL WriteVTKPoint('point-cc.vtk',m%e(p%element)%xc)
        !CALL WriteVTKPoint('point-a.vtk',p%rOld)
        !CALL WriteVTKPoint('point-b.vtk',p%r)

        !PRINT *, "Stopping in multi found"
        !CALL stopProgram(1)
    
      END IF

!-----------------------------------------------------------      
!     Get local tet con, xez and shape functions
!----------------------------------------------------------- 

      DO i = 1,e%nVertex

        IF  (e%nodecon(i).EQ.e%tet(ifound)%con(1) ) THEN
          loctcon(1) = i
        ELSE IF  (e%nodecon(i).EQ.e%tet(ifound)%con(2) ) THEN
          loctcon(2) = i
        ELSE IF  (e%nodecon(i).EQ.e%tet(ifound)%con(3) ) THEN
          loctcon(3) = i
        END IF 

      END DO

      x1 = m%x( e%tet(ifound)%con(1),: )
      x2 = m%x( e%tet(ifound)%con(2),: )
      x3 = m%x( e%tet(ifound)%con(3),: )
      x4 = e%xc

      CALL tet_kks2lks(r,x1,x2,x3,x4,xez,dummierr)

      CALL tet_shapef(xez,fi)

!-----------------------------------------------------------      
!     Interpolate in chosen tet with shape functions
!----------------------------------------------------------- 

      !PRINT *, "nvert=",e%nVertex
      !PRINT *, "con=",e%tet(ifound)%con
      !PRINT *, "lcon=",loctcon
      !PRINT *, "fi=",fi

      !
      !interpolate :: v
      !
      DO i = 1,e%nVertex
        inode = e%nodeCon(i)
        uin1(i) = fluid%Un(inode,1)
        uin2(i) = fluid%Un(inode,2)
        uin3(i) = fluid%Un(inode,3)
      END DO

      fap%vx = interpolate(e%nVertex, loctcon, uin1, fi)
      fap%vy = interpolate(e%nVertex, loctcon, uin2, fi)
      fap%vz = interpolate(e%nVertex, loctcon, uin3, fi)

      !
      !interpolate :: w
      !
      DO i = 1,e%nVertex
        inode = e%nodeCon(i)
        uin1(i) = fluid%Vortn(inode,1)
        uin2(i) = fluid%Vortn(inode,2)
        uin3(i) = fluid%Vortn(inode,3)
      END DO

      fap%wx = interpolate(e%nVertex, loctcon, uin1, fi)
      fap%wy = interpolate(e%nVertex, loctcon, uin2, fi)
      fap%wz = interpolate(e%nVertex, loctcon, uin3, fi)

      !
      !interpolate :: gradV
      !
      DO j = 1,3

        DO i = 1,e%nVertex
          inode = e%nodeCon(i)
          uin1(i) = fluid%gradUxn(inode,j)
          uin2(i) = fluid%gradUyn(inode,j)
          uin3(i) = fluid%gradUzn(inode,j)
        END DO

        fap%gradVx(j) = interpolate(e%nVertex, loctcon, uin1, fi)
        fap%gradVy(j) = interpolate(e%nVertex, loctcon, uin2, fi)
        fap%gradVz(j) = interpolate(e%nVertex, loctcon, uin3, fi)

      END DO

      !
      !interpolate :: dvdt
      !
      DO i = 1,e%nVertex
        inode = e%nodeCon(i)
        uin1(i) = fluid%dvxdt(inode)
        uin2(i) = fluid%dvydt(inode)
        uin3(i) = fluid%dvzdt(inode)
      END DO

      fap%dvxdt = interpolate(e%nVertex, loctcon, uin1, fi)
      fap%dvydt = interpolate(e%nVertex, loctcon, uin2, fi)
      fap%dvzdt = interpolate(e%nVertex, loctcon, uin3, fi)

12    CONTINUE

      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      FUNCTION interpolate(np,con,uin,fi) RESULT(uout)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: np
        INTEGER, INTENT(IN) :: con(3)
        REAL(8), INTENT(IN) :: uin(np)
        REAL(8), INTENT(IN) :: fi(4)
        
        INTEGER i,inode
        REAL(8) ucent
        REAL(8) uout

        !get element center value
        ucent = 0.0D0
        DO i = 1,np
          ucent = ucent + uin(i)
        END DO

        ucent = ucent / np

        uout = 0.0D0
        DO i = 1,3
          inode = con(i)
          uout = uout + fi(i) * uin(inode)
        END DO

        !i=4 - center value
        uout = uout + fi(4) * ucent


    
      END FUNCTION interpolate

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE GetXiEtaZeta(m,e,r,xez,ierr)
!
!     Find local coordinates from (x,y,z) 
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      USE logFile      
      IMPLICIT NONE 

      TYPE(elementType) :: e  ! mesh element
      TYPE(meshType) :: m  ! mesh data structure
      CHARACTER(255) vrstica
      REAL(8) xez(3),pripoints(6,3),pyrpoints(5,3)
      REAL(8) x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xi,eta,zeta
      INTEGER ierr,i,j
      REAL(8) r(3) ! point inside element
      REAL(8) p(3) ! fifth point of pyramid
      REAL(8) n(3) ! normal to square face
      REAL(8) c(3) ! center of square face      

      IF (e%type.EQ.12) THEN ! hexa
        x=r(1)
        y=r(2)
        z=r(3)

        x1=m%x(e%con(1,1),1)
        y1=m%x(e%con(1,1),2)
        z1=m%x(e%con(1,1),3)

        x2=m%x(e%con(1,2),1)
        y2=m%x(e%con(1,2),2)
        z2=m%x(e%con(1,2),3)

        x4=m%x(e%con(1,3),1)
        y4=m%x(e%con(1,3),2)
        z4=m%x(e%con(1,3),3)

        x3=m%x(e%con(1,4),1)
        y3=m%x(e%con(1,4),2)
        z3=m%x(e%con(1,4),3)

        x5=m%x(e%con(1,5),1)
        y5=m%x(e%con(1,5),2)
        z5=m%x(e%con(1,5),3)

        x6=m%x(e%con(1,6),1)
        y6=m%x(e%con(1,6),2)
        z6=m%x(e%con(1,6),3)

        x8=m%x(e%con(1,7),1)
        y8=m%x(e%con(1,7),2)
        z8=m%x(e%con(1,7),3)

        x7=m%x(e%con(1,8),1)
        y7=m%x(e%con(1,8),2)
        z7=m%x(e%con(1,8),3)

        CALL hex_kks2lks(x,y,z, &
     &  x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xi,eta,zeta,ierr)

        xez(1)=xi 
        xez(2)=eta
        xez(3)=zeta

        IF (ierr.NE.0) PRINT *, "HEX,xi-eta-zeta ierr=1"

      ELSE IF (e%type.EQ.10) THEN ! tet
         CALL tet_kks2lks(r,m%x(e%con(1,1),:),m%x(e%con(1,2),:),m%x(e%con(1,3),:),m%x(e%con(1,4),:),xez,ierr)
         IF (ierr.NE.0) PRINT *, "TET,xi-eta-zeta ierr=1"
      ELSE IF (e%type.EQ.14) THEN ! pyr
        DO i=1,5
          DO j=1,3
            pyrpoints(i,j)=m%x(e%con(1,i),j)
          END DO
        END DO
        xez=0.0D0
        CALL pyr_kks2lks(r,pyrpoints,xez,ierr)  ! find xi,eta,zeta 
        IF (ierr.NE.0) PRINT *, "PYR,xi-eta-zeta ierr=1"

      ELSE IF (e%type.EQ.13) THEN ! prism
        DO i=1,6
          DO j=1,3
            pripoints(i,j)=m%x(e%con(1,i),j)
          END DO
        END DO
        CALL pri_kks2lks(r,pripoints,xez,ierr)  
        IF (ierr.NE.0) PRINT *, "PRI,xi-eta-zeta ierr=1"      
      ELSE
        WRITE (vrstica,'(A,I0,A)') "ERROR :: GetXiEtaZeta :: Element type not supported: ",e%type,"!"
        CALL logWrite (vrstica)
        CALL StopProgram(1)       
      END IF


!      if (ierr.NE.0) THEN
!            WRITE (vrstica,'(A,I0,A)') "WARNING :: GetXiEtaZeta :: xi,eta,zeta not found in el. type: ",e%type,"!"
!            CALL logWrite (vrstica)
!            print *,"centers"
!            DO i=1,e%nSide
!                  print *,i,e%centers(i,1),e%centers(i,2),e%centers(i,3)
!            END DO
!            print *,"r===",r
!            DO i=1,3
!                  IF (xez(i).LT.-1.0D0) xez(i)=-1.0D0
!                  IF (xez(i).GT. 1.0D0) xez(i)= 1.0D0
!            END DO
!            ierr=0
!      END IF

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE FindParticleInUnstructuredMesh1(m,p,iT,ierr)
!
!     Finds element in mesh, where the particle is located
!     Ref: Macpherson, Nordin and Weller. 2008. Particle Tracking in unstructured,
!     arbitrary polyhedral meshes for use in CFD and molecular dynamics          
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) :: p ! particle
      TYPE(meshType) :: m  ! mesh data structure      
      INTEGER ierr,iT,ie,ieold

      INTEGER nf !number of faces
      INTEGER nhit !number of faces hit
      REAL(8) a(3),b(3)

      INTEGER k,kmax

      REAL(8), ALLOCATABLE :: lambdaCs(:,:)
      REAL(8), ALLOCATABLE :: lambdaAs(:,:)

      INTEGER, ALLOCATABLE :: faceCon(:)
      
      INTEGER imin,iface
      REAL(8) lambdaA

      ierr = 1
      kmax = 1000  ! break search after kmax
      ieold = 0

!
!     First find lambdaC for all faces
!      
      ie = p%element
      k=0

11    CONTINUE

      k=k+1

      IF (k.GT.kmax) GOTO 12
      IF (ie.LE.0) GOTO 12

      ieold = ie

      IF ( ALLOCATED(faceCon) ) DEALLOCATE( faceCon )
      IF ( ALLOCATED(lambdaCs) ) DEALLOCATE( lambdaCs )
      IF ( ALLOCATED(lambdaAs) ) DEALLOCATE( lambdaAs )

      nf = size( m%e(ie)%faceCon )

      ALLOCATE( lambdaCs(nf,2) )
      ALLOCATE( faceCon(nf) )
      
      faceCon = m%e(ie)%faceCon ! pass element ie face con
      a = m%e(ie)%xc ! a= element center point, for lambdaC
      b = p%r
 
      CALL PartTrackFindLambdas(m,a,b,faceCon,lambdaCs,nf,nhit)
      IF (nhit.GT.0 ) THEN
        PRINT *,"--------------------------"
        WRITE(*,'(A10,I6,A4,I6)') "iT=",iT,"k=",k
        WRITE(*,'(A10,I6)') "ie=",ie
        WRITE(*,'(A10,3F12.6)') "ac=",a
        WRITE(*,'(A10,3F12.6)') "b=",b
        WRITE(*,'(A10,100I6)') "fcon",faceCon
        WRITE(*,'(A10,100F12.4)') "lmbd=",lambdaCs
        WRITE(*,'(A10,I6)') "nf=",nf
        WRITE(*,'(A10,I6)') "nhit=",nhit
      END IF

!
!     Proceed only if we have a hit, else particle stays in present ie
!
      IF ( nhit.GT.0 ) THEN

        DEALLOCATE( faceCon )
        ALLOCATE( faceCon(nhit) )
        ALLOCATE( lambdaAs(nhit,2) )

        !Find lambdaAs for entire new faceCon, passed from lambdaC
        faceCon = INT( lambdaCs(1:nhit,2) )    
        a = p%rOld
        b = p%r
        nf = nhit

        CALL PartTrackFindLambdas(m,a,b,faceCon,lambdaAs,nf,nhit) 
        PRINT *,"-----lambdaA now------"
        WRITE(*,'(A10,3F12.6)') "a=",a
        WRITE(*,'(A10,3F12.6)') "b=",b
        WRITE(*,'(A10,I6)') "nf=",nf
        WRITE(*,'(A10,I6)') "nhit=",nhit
        WRITE(*,'(A10,10I6)') "fcon2",faceCon
!
!       If we have a hit, find min lambdaA, and set new element
!
        IF (nhit.GT.0 ) THEN
          PRINT *, "yeshit:",nhit
          imin = MINLOC( lambdaAs(:,1), DIM=1, MASK=(lambdaAs(:,1) > 0) )

          lambdaA = lambdaAs(imin,1)
          iface = INT (lambdaAs(imin,2))

          WRITE(*,'(A10,F12.4)') "lmbdA=",lambdaA
          WRITE(*,'(A10,I6)') "ifaceA=",iface

          ie = setNewIe(iface)

          !WRITE(*,'(A10,3F14.6)') "aold=",a

          !move a to crosssection point
          a = a + lambdaA * ( b - a )

          !WRITE(*,'(A10,3F14.6)') "anew=",a

          GOTO 11

        ELSE
          PRINT *, "nohit:",nhit
          WRITE(*,'(A10,10F12.4)') "lmbdAs=",lambdaAs

          IF (nf.EQ.1)  THEN
            iface = INT (lambdaCs(1,2))
            ie = setNewIe(iface)

          ELSE !cant be zero, so if here then nf > 1
           
!           More than 1 face identified by Cell center test,
!           but none of the faces is crossed by |ab| (lambdaA is not in range 0...1)
!
            PRINT *, "Opet problem, nf .NE. 1"
            !CALL stopProgram(1)

          ! Return ierr = 1 for now, I don't know what to do yet, 
          ! this is getting really annoying... 

            !ierr = 1
            !RETURN

          END IF
          
        END IF


      END IF


      IF ( ALLOCATED(faceCon) ) DEALLOCATE( faceCon )
      IF ( ALLOCATED(lambdaCs) ) DEALLOCATE( lambdaCs )
      IF ( ALLOCATED(lambdaAs) ) DEALLOCATE( lambdaAs )

      p%element = ie
      ierr = 0

      RETURN

12    CONTINUE

      PRINT *, "At FindPartInUnstr.Mesh, GOTO 12"
      PRINT *, "ie=",ie
      PRINT *, "k=",k,"kmax=",kmax

      IF (k.GT.kmax) CALL stopProgram(1)


      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      !
      ! Finds which index in array contains lambdaMin
      !

      INTEGER FUNCTION setNewIe(iface)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: iface

        INTEGER ie

        !PRINT *, "ifce=",iface
        !PRINT *, "own=",m%f(abs(iface))%owner
        !PRINT *, "nei=",m%f(abs(iface))%neighbour

        IF (iface.GT.0) THEN !old element was face owner
          ie = m%f(iface)%neighbour
        ELSE IF (iface.LT.0) THEN !old element was face neighobur
          ie = m%f(-iface)%owner
        ELSE
          PRINT *, "Error :: iface = 0"
          CALL stopProgram(1)
        END IF 

        setNewIe = ie
    
      END FUNCTION setNewIe


      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE FindParticleInUnstructuredMesh2(m,p,iT,ierr)
!
!     Finds element in mesh, where the particle is located
!     Ref: Macpherson, Nordin and Weller. 2008. Particle Tracking in unstructured,
!     arbitrary polyhedral meshes for use in CFD and molecular dynamics          
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) :: p ! particle
      TYPE(meshType) :: m  ! mesh data structure      
      INTEGER ierr,iT,ie,ieold

      INTEGER nf !number of faces
      INTEGER nhit !number of faces hit
      REAL(8) a(3),b(3)

      INTEGER i,j,k,kmax

      INTEGER iehist(5)

      REAL(8), ALLOCATABLE :: lambdaCs(:)
      REAL(8), ALLOCATABLE :: lambdaAs(:)

      INTEGER, ALLOCATABLE :: faceCon(:)
      INTEGER, ALLOCATABLE :: faceCon2(:)
      INTEGER, ALLOCATABLE :: hitArr(:)
      
      INTEGER imin,iface
      REAL(8) lambdaA,lambdaM

      ierr = 1
      kmax = 4  ! break search after kmax
      ieold = 0
      iehist = 0

!
!     First find lambdaC for all faces
!      
      ie = p%element
      k=0

11    CONTINUE

      k=k+1

      !Crossed boundary face
      IF (ie.LE.0) GOTO 12
      
      !Big error
      IF (k.GT.kmax) GOTO 12

      !Proceed without detection
      DO i=1,k
        IF (ie.EQ.iehist(i)) THEN
          PRINT *, "Error :: same ie twice"
          ierr = 2
          RETURN
        END IF
      END DO

      iehist(k) = ie
      ieold = ie

      IF ( ALLOCATED(faceCon) ) DEALLOCATE( faceCon )
      IF ( ALLOCATED(faceCon2) ) DEALLOCATE( faceCon2 )
      IF ( ALLOCATED(hitArr) ) DEALLOCATE( hitArr )
      IF ( ALLOCATED(lambdaCs) ) DEALLOCATE( lambdaCs )
      IF ( ALLOCATED(lambdaAs) ) DEALLOCATE( lambdaAs )

      nf = size( m%e(ie)%faceCon )

      ALLOCATE( lambdaCs(nf) )
      ALLOCATE( faceCon(nf) )
      ALLOCATE( hitArr(nf) )
      
      faceCon = m%e(ie)%faceCon ! pass element ie face con
      a = m%e(ie)%xc ! a= element center point, for lambdaC
      b = p%r
 
      CALL PartTrackFindLambdas2(m,a,b,faceCon,lambdaCs,hitArr,nf,nhit)
      IF (nhit.GT.0 ) THEN
        PRINT *,"--------------------------"
        WRITE(*,'(A10,I6,A4,I6)') "iT=",iT,"k=",k
        WRITE(*,'(A10,I6)') "ie=",ie
        WRITE(*,'(A10,3F12.6)') "ac=",a
        WRITE(*,'(A10,3F12.6)') "b=",b
        WRITE(*,'(A10,100I6)') "fcon",faceCon
        WRITE(*,'(A10,100F12.4)') "lmbd=",lambdaCs
        WRITE(*,'(A10,I6)') "nf=",nf
        WRITE(*,'(A10,I6)') "nhit=",nhit
      END IF

!
!     Proceed only if we have a hit, else particle stays in present ie
!
      IF ( nhit.GT.0 ) THEN

        ALLOCATE( faceCon2(nhit) )
      
        j=0
        DO i = 1,nf
            IF ( hitArr(i).EQ.1 ) THEN
              j = j + 1
              faceCon2(j)=faceCon(i)
            END IF
        END DO

        DEALLOCATE( faceCon )
        DEALLOCATE( hitArr )
        ALLOCATE( hitArr(nhit) )
        ALLOCATE( lambdaAs(nhit) )

        !Find lambdaAs for entire new faceCon, passed from lambdaC  
        a = p%rOld
        b = p%r
        nf = nhit

        CALL PartTrackFindLambdas2(m,a,b,faceCon2,lambdaAs,hitArr,nf,nhit) 
        PRINT *,"-----lambdaA now------"
        WRITE(*,'(A10,3F12.6)') "a=",a
        WRITE(*,'(A10,3F12.6)') "b=",b
        WRITE(*,'(A10,I6)') "nf=",nf
        WRITE(*,'(A10,I6)') "nhit=",nhit
        WRITE(*,'(A10,10I6)') "fcon2",faceCon

        !Find min lambdaA, and change ie if neccessary
        imin = MINLOC(lambdaAs,DIM=1)
        lambdaA = lambdaAs(imin)

        !do ref eq. 5: p = a+lm*(b-a); lm = min(1,max(0,la))
        lambdaM = MINVAL((/1.0D0, MAXVAL((/0.0D0,lambdaA/)) /))

        !move a to crosssection point
        WRITE(*,'(A10,3F14.6)') "aold=",a
        a = a + lambdaM * ( b - a )

        !set cell occupany according to definitions in ref article
        iface = faceCon2(imin)
        ie = setNewIe(iface)

        WRITE(*,'(A10,3F12.4)') "lmbdAs=",lambdaAs
        WRITE(*,'(A10,F12.4)') "lmbdM=",lambdaM
        WRITE(*,'(A10,I6)') "ifaceA=",iface        
        WRITE(*,'(A10,3F14.6)') "anew=",a
        WRITE(*,'(A10,I6)') "ienew=",ie

        !go back in loop only for lambda < 1
        IF (lambdaA.LT.0) GOTO 11

      END IF


      IF ( ALLOCATED(faceCon) ) DEALLOCATE( faceCon )
      IF ( ALLOCATED(faceCon2) ) DEALLOCATE( faceCon2 )
      IF ( ALLOCATED(hitArr) ) DEALLOCATE( hitArr )
      IF ( ALLOCATED(lambdaCs) ) DEALLOCATE( lambdaCs )
      IF ( ALLOCATED(lambdaAs) ) DEALLOCATE( lambdaAs )


      !Check again for bad ie
      IF (ie.LE.0) GOTO 12

    
      p%element = ie

13    CONTINUE 
      ierr = 0

      RETURN

12    CONTINUE

      PRINT *, "At FindPartInUnstr.Mesh, GOTO 12"
      PRINT *, "ie=",ie
      PRINT *, "k=",k,"kmax=",kmax

      !IF (k.GT.kmax) CALL stopProgram(1)


      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      !
      ! Finds which index in array contains lambdaMin
      !

      INTEGER FUNCTION setNewIe(iface)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: iface

        INTEGER ie

        !PRINT *, "ifce=",iface
        !PRINT *, "own=",m%f(abs(iface))%owner
        !PRINT *, "nei=",m%f(abs(iface))%neighbour

        IF (iface.GT.0) THEN !old element was face owner
          ie = m%f(iface)%neighbour
        ELSE IF (iface.LT.0) THEN !old element was face neighobur
          ie = m%f(-iface)%owner
        ELSE
          PRINT *, "Error :: iface = 0"
          CALL stopProgram(1)
        END IF 

        setNewIe = ie
    
      END FUNCTION setNewIe


      END      

! -----------------------------------------------------------------------------------------
      SUBROUTINE FindParticleInUnstructuredMesh3(m,p,iT,ierr)
!
!     Finds element in mesh, where the particle is located
!     Ref: Macpherson, Nordin and Weller. 2008. Particle Tracking in unstructured,
!     arbitrary polyhedral meshes for use in CFD and molecular dynamics          
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) :: p ! particle
      TYPE(meshType) :: m  ! mesh data structure      
      INTEGER ierr,iT,ie

      INTEGER nf,nt
      INTEGER nhit !number hits

      REAL(8) a(3),b(3),x(3,3),fb(3)
    
      INTEGER i,j,k,l
      INTEGER iface
      REAL(8) cc(3),fc(3),fn(3)

      INTEGER, ALLOCATABLE :: fcon(:,:)
      REAL(8), ALLOCATABLE :: centers(:,:), normals(:,:),lambdas(:,:)


      INTEGER, ALLOCATABLE :: fcon2(:,:)
      REAL(8), ALLOCATABLE :: centers2(:,:), normals2(:,:),lambdas2(:,:)

      INTEGER :: lastHit(2)

      lastHit = 0
 
      ie = p%element 
      
      a = p%rOld
      b = p%r

10    CONTINUE

      IF (ie.LE.0) GOTO 11

      IF ( ALLOCATED(fcon) ) DEALLOCATE( fcon )
      IF ( ALLOCATED(centers) ) DEALLOCATE( centers )
      IF ( ALLOCATED(normals) ) DEALLOCATE( normals )
      IF ( ALLOCATED(lambdas) ) DEALLOCATE( lambdas )

      IF ( ALLOCATED(fcon2) ) DEALLOCATE( fcon2 )
      IF ( ALLOCATED(centers2) ) DEALLOCATE( centers2 )
      IF ( ALLOCATED(normals2) ) DEALLOCATE( normals2 )
      IF ( ALLOCATED(lambdas2) ) DEALLOCATE( lambdas2 )

      !Count number of all triangles in element surface - ntall
      nf = size ( m%e(ie)%faceCon )

      k=0
      DO i = 1,nf

        iface = ABS( m%e(ie)%faceCon(i) )
        nt = size( m%f(iface)%tf )
        k = k + nt

      END DO

      ALLOCATE( fcon(k,2) )
      ALLOCATE( centers(k,3) )
      ALLOCATE( normals(k,3) )
      ALLOCATE( lambdas(k,2) )

      !Prepare connectivity, centers and normals
      k=0
      DO i = 1,nf

        iface = abs( m%e(ie)%faceCon(i) )

        nt = size( m%f(iface)%tf )
        DO j = 1,nt

          k = k+1
          fcon(k,1) = m%e(ie)%faceCon(i)
          fcon(k,2) = j

          centers(k,:) = m%f(iface)%tf(j)%center
          normals(k,:) = m%f(iface)%tf(j)%normal

        END DO

      END DO

      !IF (nhit.EQ.0) GOTO 11
!
      !ALLOCATE( fcon2(nhit,2) )
      !ALLOCATE( centers2(nhit,3) )
      !ALLOCATE( normals2(nhit,3) )
      !ALLOCATE( lambdas2(nhit,2) )
!
      !DO i = 1,k
      !  IF (lambdas(i,2).GT.0.5D0) THEN
      !    fcon2(i,:) = fcon(i,:)
      !    centers2(i,:) = centers(i,:)
      !    normals2(i,:) = normals(i,:)
      !  END IF
      !END DO

      !Get face lambdas
      !PRINT *, "getting lambda for p%id=",p%id
      CALL PartTrackFindLambdas3(a,b,centers,normals,k,lambdas,nhit)

      IF (nhit.GT.0) THEN

        l = 0
        DO i = 1,k 

          !WRITE(*,'(A5,2I6,A4,F16.4)') "fcon=",fcon(i,:)," L=",lambdas(i,1)

          IF (lambdas(i,2).GT.0.5D0) THEN

            a = a + lambdas(i,1) * ( b - a )

            !Test barycentric coordinates

            iface = ABS( fcon(i,1) )
            j = fcon(i,2)

            !Skip if this is the last face and triangle hit
            !Happens due to numerical errors
            IF(iface.EQ.lastHit(1).AND.j.EQ.lastHit(2)) EXIT

            x(1,:) = m%x( m%f(iface)%tf(j)%con(1),: )
            x(2,:) = m%x( m%f(iface)%tf(j)%con(2),: )
            x(3,:) = m%x( m%f(iface)%tf(j)%con(3),: )

            fc = m%f(iface)%tf(j)%center
            fn = m%f(iface)%tf(j)%normal

            CALL GetBarycentricInTri(a,x,fc,fn,fb,ierr)
            WRITE(*,'(A4,3F12.6)') "fb=",fb
            PRINT *, "ierr=",ierr

            IF (ierr.EQ.0) THEN

              l = l + 1 
              ie = setNewIe( fcon(i,1) )

              lastHit(1) = abs(fcon(i,1))
              lastHit(2) = fcon(i,2)

            END IF

          END IF

        END DO

        

        IF (l.EQ.1) THEN
          GOTO 10
        ELSE IF (l.GT.1) THEN
          PRINT *, "More than one hit, l=",l 
          CALL StopProgram(1) 

        END IF

      END IF

11    CONTINUE


      !Return element
      IF (ie.GT.0) THEN

        p%element = ie
        ierr = 0

      !Crossed boundary face or other problem
      ELSE 

        ierr = 1

      END IF

!-------------------------------------------------------------


!      IF (ie.LE.0) GOTO 11
!
!      IF ( ALLOCATED(fcon) ) DEALLOCATE( fcon )
!      IF ( ALLOCATED(centers) ) DEALLOCATE( centers )
!      IF ( ALLOCATED(normals) ) DEALLOCATE( normals )
!      IF ( ALLOCATED(lambdas) ) DEALLOCATE( lambdas )
!
!      IF ( ALLOCATED(fcon2) ) DEALLOCATE( fcon2 )
!      IF ( ALLOCATED(centers2) ) DEALLOCATE( centers2 )
!      IF ( ALLOCATED(normals2) ) DEALLOCATE( normals2 )
!      IF ( ALLOCATED(lambdas2) ) DEALLOCATE( lambdas2 )
!
!      !Count number of all triangles in element surface - ntall
!      nf = size ( m%e(ie)%faceCon )
!
!      k=0
!      DO i = 1,nf
!
!        iface = ABS( m%e(ie)%faceCon(i) )
!        nt = size( m%f(iface)%tf )
!        k = k + nt
!
!      END DO
!
!      ALLOCATE( fcon(k,2) )
!      ALLOCATE( centers(k,3) )
!      ALLOCATE( normals(k,3) )
!      ALLOCATE( lambdas(k,2) )
!
!      !Prepare connectivity, centers and normals
!      k=0
!      DO i = 1,nf
!
!        iface = abs( m%e(ie)%faceCon(i) )
!
!        nt = size( m%f(iface)%tf )
!        DO j = 1,nt
!
!          k = k+1
!          fcon(k,1) = m%e(ie)%faceCon(i)
!          fcon(k,2) = j
!
!          centers(k,:) = m%f(iface)%tf(j)%center
!          normals(k,:) = m%f(iface)%tf(j)%normal
!
!        END DO
!
!      END DO
!
!      !!Get lambdaC's 
!      cc = m%e(ie)%xc
!      CALL PartTrackFindLambdas3(cc,b,centers,normals,k,lambdas,nhit)
!
!      IF (nhit.GT.0) THEN
!        PRINT *, "LambdaC nhit=",nhit
!        PRINT *, "p%element=",p%element 
!        PRINT *, "ie=",ie
!        PRINT *, "Stopping program"
!        CALL StopProgram(1)
!      END IF





      !-------------------------------------------------------------



      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      !
      ! Finds which index in array contains lambdaMin
      !

      INTEGER FUNCTION setNewIe(iface)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: iface

        INTEGER ie

        !PRINT *, "ifce=",iface
        !PRINT *, "own=",m%f(abs(iface))%owner
        !PRINT *, "nei=",m%f(abs(iface))%neighbour

        IF (iface.GT.0) THEN !old element was face owner
          ie = m%f(iface)%neighbour
        ELSE IF (iface.LT.0) THEN !old element was face neighobur
          ie = m%f(-iface)%owner
        ELSE
          PRINT *, "Error :: iface = 0"
          CALL stopProgram(1)
        END IF 

        setNewIe = ie
    
      END FUNCTION setNewIe


      END 


! -----------------------------------------------------------------------------------------
      SUBROUTINE FindParticleInUnstructuredMesh4(m,p,iT,ierr,debugout)
!
!     Finds element in mesh, where the particle is located
!     Ref: Macpherson, Nordin and Weller. 2008. Particle Tracking in unstructured,
!     arbitrary polyhedral meshes for use in CFD and molecular dynamics          
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) :: p ! particle
      TYPE(meshType) :: m  ! mesh data structure      
      INTEGER ierr,iT,ie

      INTEGER nf,nt
      INTEGER nhit !number hits

      REAL(8) a(3),b(3),ac(3),x(3,3),fb(3)
    
      INTEGER i,j,k,l,h
      INTEGER iface
      REAL(8) cc(3),fc(3),fn(3), small

      INTEGER, ALLOCATABLE :: fcon(:,:)
      REAL(8), ALLOCATABLE :: centers(:,:), normals(:,:),lambdas(:,:)


      INTEGER, ALLOCATABLE :: fcon2(:,:)
      REAL(8), ALLOCATABLE :: centers2(:,:), normals2(:,:),lambdas2(:,:)

      REAL(8), ALLOCATABLE :: baryHit(:,:)

      INTEGER :: lastHit(2),dcnt

      REAL(8) :: debugout(100,3)

      LOGICAL doLambdaC

      doLambdaC = .FALSE.
      small = 1.0E-20

      PRINT *, "--------------------------------"
      PRINT *, "New Find"
      PRINT *, "--------------------------------"

      dcnt = 10

      lastHit = 0
      h = 0
      ie = p%element 
      a = p%rOld
      b = p%r

      debugout(1,:) = a
      debugout(2,:) = b

10    CONTINUE

      debugout(dcnt,1) = ie
      dcnt = dcnt + 1

      !PRINT *, "a=",a 
      !PRINT *, "b=",b
      PRINT *,"------------------"
      PRINT *, "Current ie=",ie
      WRITE(*,'(A4,3F16.10)') "a=",a
      WRITE(*,'(A4,3F16.10)') "b=",b

      IF (ie.LE.0) GOTO 11

      !PRINT *, "Starting particle find"

      IF ( ALLOCATED(fcon) ) DEALLOCATE( fcon )
      IF ( ALLOCATED(centers) ) DEALLOCATE( centers )
      IF ( ALLOCATED(normals) ) DEALLOCATE( normals )
      IF ( ALLOCATED(lambdas) ) DEALLOCATE( lambdas )

      IF ( ALLOCATED(fcon2) ) DEALLOCATE( fcon2 )
      IF ( ALLOCATED(centers2) ) DEALLOCATE( centers2 )
      IF ( ALLOCATED(normals2) ) DEALLOCATE( normals2 )
      IF ( ALLOCATED(lambdas2) ) DEALLOCATE( lambdas2 )

      !Count number of all triangles in element surface - ntall
      nf = size ( m%e(ie)%faceCon )

      k=0
      DO i = 1,nf

        iface = ABS( m%e(ie)%faceCon(i) )
        nt = size( m%f(iface)%tf )
        k = k + nt

      END DO

      ALLOCATE( fcon(k,2) )
      ALLOCATE( centers(k,3) )
      ALLOCATE( normals(k,3) )
      ALLOCATE( lambdas(k,2) )

      !Prepare connectivity, centers and normals
      k=0
      DO i = 1,nf

        iface = abs( m%e(ie)%faceCon(i) )

        nt = size( m%f(iface)%tf )
        DO j = 1,nt

          k = k+1
          fcon(k,1) = m%e(ie)%faceCon(i)
          fcon(k,2) = j

          centers(k,:) = m%f(iface)%tf(j)%center
          normals(k,:) = m%f(iface)%tf(j)%normal

        END DO

      END DO

      !Get face lambdas
      !PRINT *, "ie=",ie
      !PRINT *, "getting lambda for p%id=",p%id
            
      CALL PartTrackFindLambdas3(a,b,centers,normals,k,lambdas,nhit)
      PRINT *, "lambdaA nhits=",nhit

      IF (nhit.GT.0) THEN

        IF ( ALLOCATED(baryHit) ) DEALLOCATE(baryHit)
        ALLOCATE( baryHit(k,8) )
        baryHit = 0.0D0

        l=0
        DO i = 1,k 

          !WRITE(*,'(A5,2I6,A4,F16.4)') "fcon=",fcon(i,:)," L=",lambdas(i,1)

          IF (lambdas(i,2).GT.0.5D0) THEN

            !PRINT *, "lambda hit at i=",i

            !get face crossection point
            ac = a + lambdas(i,1) * ( b - a )

            !Get Barycentric coordinates
            iface = ABS( fcon(i,1) )
            j = fcon(i,2)

            !Skip if this is the last face and triangle hit
            !Happens due to numerical errors
            IF(iface.EQ.lastHit(1).AND.j.EQ.lastHit(2)) THEN 

              WRITE(*,'(A8,I4,A3,G16.8,A35,2I8)') "lambdaA(",i,")=", &
              & lambdas(i,1)," terminating, it's the same one",iface,j

              !PRINT *, "Terminating, due to last hit"
              !PRINT *, "iface=",iface," lasthit1=",lastHit(1)
              !PRINT *, "j=",j," lasthit2=",lastHit(2)

              !EXIT  --->> !!!!! Warning, this kills everything !!!!!!
              GOTO 14

            END IF

            x(1,:) = m%x( m%f(iface)%tf(j)%con(1),: )
            x(2,:) = m%x( m%f(iface)%tf(j)%con(2),: )
            x(3,:) = m%x( m%f(iface)%tf(j)%con(3),: )

            fc = m%f(iface)%tf(j)%center
            fn = m%f(iface)%tf(j)%normal

            !WRITE(*,'(A5,3F22.16)') "a=",a
            !WRITE(*,'(A5,3F22.16)') "b=",b

            !CALL GetBarycentricInTri(ac,x,fc,fn,fb,ierr)
            CALL GetBarycentricInTriDebug(ac,x,fc,fn,fb,ierr)

            !WRITE(*,'(A4,3F12.6)') "fb=",fb
            !PRINT *, "ierr=",ierr

            !PRINT *, "fcon=",fcon(:,1)
            !WRITE(*,'(A8,I4,A3,G16.8)') "lambdaA(",i,")=",lambdas(i,1)
            !WRITE(*,'(A4,3F12.6,G20.12,I3)') "fb=",fb,sum(fb),ierr

            IF (ierr.EQ.0) THEN

              !PRINT *, "Inside bary hit, with l=",l+1
              !WRITE(*,'(A8,I4,A3,G16.8)') "lambdaA(",i,")=",lambdas(i,1)
              !WRITE(*,'(A4,3F12.6,G20.12)') "fb=",fb,sum(fb)

              !PRINT *, "iebefore=",ie

              ie = setNewIe( fcon(i,1) )

              !PRINT *, "ieafter",ie

              !lastHit(1) = abs(fcon(i,1))
              !lastHit(2) = fcon(i,2)

              !PRINT *, "ie=",ie," chg lstht=",lastHit

              l = l + 1

              p%potentialElements(l) = ie  
              
              baryHit(l,1) = ie 
              baryHit(l,2) = lambdas(i,1)
              baryHit(l,3) = sum(fb)
              baryHit(l,4:6) = ac 
              baryHit(l,7) = abs(fcon(i,1)) !lasthit 1
              baryhit(l,8) = fcon(i,2) !last hit 2

            END IF

            WRITE(*,'(A8,I2,A3,G16.8,A4,3F10.4,G16.8,I3,A4,I10)') "lambdaA(",i,")=", &
            & lambdas(i,1)," fb=",fb,sum(fb)-1.0D0,ierr," ie=",ie

            WRITE(*,'(A6,9F18.10,A6,2I8)') "xfac=",x(1,:),x(2,:),x(3,:)," if,j=",iface,j



          ELSE !If no hit just debug write

            WRITE(*,'(A8,I4,A3,G16.8,A8)') "lambdaA(",i,")=",lambdas(i,1)," no hit"

          END IF

14        CONTINUE

        END DO      

        IF (l.GT.0) THEN

          !Find the max lambda in array
          !i = MAXLOC(baryHit(:,2),DIM=1)

          !Find the min lambda in array
          i = MINLOC(baryHit(:,2),DIM=1,MASK=(baryHit(:,2) > small) )

          !select the correct ie
          ie = INT( baryHit(i,1) )

          !change a point to ac (intersection with face)
          a = baryHit(i,4:6)

          !Fill last hit

          lastHit(1) = INT( baryHit(i,7) )
          lastHit(2) = INT( baryHit(i,8) )

          !debugging - record ie change
          debugout(dcnt,1) = ie
          dcnt = dcnt + 1

          PRINT *, "Selecting baryhit l=",i 
          PRINT *, "Selecting element ie=",ie
          WRITE(*,'(A6,3F20.10)') "anew=",a
          PRINT *, "Lasthit=",lasthit

          IF (l.GT.1) THEN

            PRINT *, "Problem, more than one hit, l=",l 
            PRINT *, "h is now=",h
            PRINT *, "------------"

          END IF

          IF (h.GT.20) THEN 
            PRINT *, "Problem :: h > hmax; h=",h 
            CALL BreakPoint() 
            GOTO 11
          END IF

          h = h + 1

          GOTO 10

        END IF

      END IF

11    CONTINUE

      IF (ie.GT.0 .AND. doLambdaC) THEN

        IF ( ALLOCATED(fcon) ) DEALLOCATE( fcon )
        IF ( ALLOCATED(centers) ) DEALLOCATE( centers )
        IF ( ALLOCATED(normals) ) DEALLOCATE( normals )
        IF ( ALLOCATED(lambdas) ) DEALLOCATE( lambdas )

        !Count number of all triangles in element surface - ntall
        nf = size ( m%e(ie)%faceCon )

        k=0
        DO i = 1,nf

          iface = ABS( m%e(ie)%faceCon(i) )
          nt = size( m%f(iface)%tf )
          k = k + nt

        END DO

        ALLOCATE( fcon(k,2) )
        ALLOCATE( centers(k,3) )
        ALLOCATE( normals(k,3) )
        ALLOCATE( lambdas(k,2) )

        !Prepare connectivity, centers and normals
        k=0
        DO i = 1,nf

          iface = abs( m%e(ie)%faceCon(i) )

          nt = size( m%f(iface)%tf )
          DO j = 1,nt

            k = k+1
            fcon(k,1) = m%e(ie)%faceCon(i)
            fcon(k,2) = j

            centers(k,:) = m%f(iface)%tf(j)%center
            normals(k,:) = m%f(iface)%tf(j)%normal

          END DO

        END DO

        cc = m%e(ie)%xc
        CALL PartTrackFindLambdas3(cc,b,centers,normals,k,lambdas,nhit)

        !PRINT *, "LambdaC nhits:",nhit

        IF (nhit.NE.0) THEN

          !CALL WriteVTKTriangulatedElement('trielement.vtk',m,ie)
          !CALL WriteVTKPoint('point-cc.vtk',m%e(ie)%xc)
          !CALL WriteVTKPoint('point-a.vtk',p%rOld)
          !CALL WriteVTKPoint('point-b.vtk',p%r)

          DO i = 1,k 

            IF (lambdas(i,2).GT.0.5D0) THEN
  
              !get face crossection point
              ac = cc + lambdas(i,1) * ( b - cc )
  
              !Get Barycentric coordinates
              iface = ABS( fcon(i,1) )
              j = fcon(i,2)
  
              x(1,:) = m%x( m%f(iface)%tf(j)%con(1),: )
              x(2,:) = m%x( m%f(iface)%tf(j)%con(2),: )
              x(3,:) = m%x( m%f(iface)%tf(j)%con(3),: )
  
              fc = m%f(iface)%tf(j)%center
              fn = m%f(iface)%tf(j)%normal
  
              CALL GetBarycentricInTri(ac,x,fc,fn,fb,ierr)
              PRINT *, "LambdaC ------"
              WRITE(*,'(A4,3F12.6,G16.8)') "fb=",fb,sum(fb)
              WRITE(*,'(A8,I4,A3,G16.8)') "lambdaC(",i,")=",lambdas(i,1)

              IF (ierr.EQ.0) THEN
                PRINT *, "LambdaC ierr = 0"
                WRITE(*,'(A4,I6)') "ie1=",ie
                ie = setNewIe( fcon(i,1) )
                WRITE(*,'(A4,I6)') "ie2=",ie
                PRINT *, "----------------------"

                !debugging - record ie change
                debugout(dcnt,1) = ie
                dcnt = dcnt + 1

                !IF (l.GT.1) CALL BreakPoint()
              END IF
  
            END IF
          END DO

        END IF

        !Set element 
        p%element = ie
        ierr = 0

      ELSE

        !Just set element, dont do lambdaC
        p%element = ie
        ierr = 0

      END IF 

      !Crossed boundary face
      IF (ie.LT.1) ierr = 2



      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      !
      ! Finds which index in array contains lambdaMin
      !

      INTEGER FUNCTION setNewIe(iface)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: iface

        INTEGER ie

        !PRINT *, "ifce=",iface
        !PRINT *, "own=",m%f(abs(iface))%owner
        !PRINT *, "nei=",m%f(abs(iface))%neighbour

        IF (iface.GT.0) THEN !old element was face owner
          ie = m%f(iface)%neighbour
        ELSE IF (iface.LT.0) THEN !old element was face neighobur
          ie = m%f(-iface)%owner
        ELSE
          PRINT *, "Error :: iface = 0"
          CALL stopProgram(1)
        END IF 

        setNewIe = ie
    
      END FUNCTION setNewIe


      END





! -----------------------------------------------------------------------------------------
      SUBROUTINE PartTrackFindLambdas(m,a,b,fcon,lambdaArr,nf,nhit)
!
!     Calculates the ratio of tracking step ("track") from start to the 
!     intersection with a face (its plane actually)
!        
!     Inputs:
!     m - mesh TYPE
!     p - particle TYPE
!     ie - element index 
!     a - particle r old (x,y,z)
!     b - particle r (x,y,z)
!        
!     Returns:
!     lambdaArr(:,1) - resulting ratios; valid = 0..1; invalid (no intersection) = <0 or >1
!     lambdaArr(:,2) - signed iface, with which ratio was calculated
!     nhit -> number of faces hit
!            
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      IMPLICIT NONE

      TYPE(meshType) :: m 

      INTEGER i,nf,iface
      REAL(8) Cf(3),a(3),b(3),n(3)
      INTEGER fcon(nf)
      REAL(8) fc(nf,3)
      REAL(8) lambda,lambdaArr(nf,2)

      INTEGER nhit !number of faces hit

      nhit = 0
      lambdaArr = 0.0D0

      DO i = 1,nf

        iface = ABS( fcon(i) )

        Cf = m%f(iface)%center
        n = m%f(iface)%normal

        !WRITE(*,'(A6,3F12.4)') " a=",a
        !WRITE(*,'(A6,3F12.4)') " b=",b 
        !WRITE(*,'(A6,3F12.4)') " n=",n
        !WRITE(*,'(A6,3F12.4)') " Cf=",Cf

        lambda = calcLambda(Cf,a,b,n)
        lambdaArr(i,1)=lambda

        IF (lambda.GT.0.0D0.AND.lambda.LE.1.0D0) THEN

          nhit = nhit+1

          lambdaArr(nhit,1)=lambda
          lambdaArr(nhit,2)=fcon(i)

        END IF

      END DO



      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------
      REAL(8) FUNCTION calcLambda(Cf,a,b,n)

        !IMPLICIT NONE
        REAL(8), INTENT(IN) :: Cf(3),a(3),b(3),n(3)

        REAL (8) lmbd, dp_up, dp_dwn

        lmbd = -1.0D0

        !numerator
        CALL DotProduct(Cf-a,n,dp_up)

        !denominator
        CALL DotProduct(b-a,n,dp_dwn)

        !if denominator not zero
        IF ( ABS(dp_dwn).GT.1.0E-20 ) THEN
            lmbd = dp_up / dp_dwn
        END IF

        calcLambda = lmbd
      
      END FUNCTION calcLambda

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE PartTrackFindLambdas2(m,a,b,fcon,lambdaArr,hitArr,nf,nhit)
!
!     Calculates the ratio of tracking step ("track") from start to the 
!     intersection with a face (its plane actually)
!        
!     Inputs:
!     m - mesh TYPE
!     p - particle TYPE
!     ie - element index 
!     a - particle r old (x,y,z)
!     b - particle r (x,y,z)
!        
!     Returns:
!            
! -----------------------------------------------------------------------------------------
      USE mesh
      USE superE
      IMPLICIT NONE

      TYPE(meshType) :: m 

      INTEGER i,nf,iface
      REAL(8) Cf(3),a(3),b(3),n(3)
      INTEGER fcon(nf)
      REAL(8) lambda

      REAL(8) lambdaArr(nf)
      INTEGER hitArr(nf)

      INTEGER nhit !number of faces hit

      nhit = 0
      lambdaArr = 0.0D0
      hitArr = 0

      DO i = 1,nf

        iface = ABS( fcon(i) )

        Cf = m%f(iface)%center
        n = m%f(iface)%normal

        !WRITE(*,'(A6,3F12.4)') " a=",a
        !WRITE(*,'(A6,3F12.4)') " b=",b 
        !WRITE(*,'(A6,3F12.4)') " n=",n
        !WRITE(*,'(A6,3F12.4)') " Cf=",Cf

        lambda = calcLambda(Cf,a,b,n)
        lambdaArr(i)=lambda

        IF (lambda.GT.0.0D0.AND.lambda.LE.1.0D0) THEN

          nhit = nhit+1
          hitArr(i) = 1

        END IF

      END DO


      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------
      REAL(8) FUNCTION calcLambda(Cf,a,b,n)

        !IMPLICIT NONE
        REAL(8), INTENT(IN) :: Cf(3),a(3),b(3),n(3)

        REAL (8) lmbd, dp_up, dp_dwn

        lmbd = -1.0D0

        !numerator
        CALL DotProduct(Cf-a,n,dp_up)

        !denominator
        CALL DotProduct(b-a,n,dp_dwn)

        !if denominator not zero
        IF ( ABS(dp_dwn).GT.1.0E-20 ) THEN
            lmbd = dp_up / dp_dwn
        END IF

        calcLambda = lmbd
      
      END FUNCTION calcLambda

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE PartTrackFindLambdas3(a,b,centers,normals,nf,lambdas,hit)
!
!     Calculates the ratio of tracking step ("track") from start to the 
!     intersection with a face (its plane actually)
!        
!     Inputs:
!     a - particle r old (x,y,z)
!     b - particle r (x,y,z)
!        
!     Returns:
!            
! -----------------------------------------------------------------------------------------

      IMPLICIT NONE

      REAL(8) a(3),b(3)
      INTEGER nf, hit
      REAL(8) centers(nf,3)
      REAL(8) normals(nf,3)
      REAL(8) lambdas(nf,2)

      REAL(8) Cf(3),n(3)

      INTEGER i

      lambdas = 0.0D0
      hit = 0

      DO i = 1,nf

        Cf = centers(i,:)
        n = normals(i,:)

        !WRITE(*,'(A6,3F12.4)') " a=",a
        !WRITE(*,'(A6,3F12.4)') " b=",b 
        !WRITE(*,'(A6,3F12.4)') " n=",n
        !WRITE(*,'(A6,3F12.4)') " Cf=",Cf

        lambdas(i,1) = calcLambda(Cf,a,b,n)

        IF(lambdas(i,1).GT.0.0D0.AND.lambdas(i,1).LE.1.0D0) THEN
          lambdas(i,2) = 1
          hit = hit + 1
        END IF

        !WRITE(*,'(A6,F12.5)') " L=",lambdas(i)

      END DO


      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------
      REAL(8) FUNCTION calcLambda(Cf,a,b,n)

        !IMPLICIT NONE
        REAL(8), INTENT(IN) :: Cf(3),a(3),b(3),n(3)

        REAL (8) lmbd, dp_up, dp_dwn

        lmbd = -1.0D0

        !numerator
        CALL DotProduct(Cf-a,n,dp_up)

        !denominator
        CALL DotProduct(b-a,n,dp_dwn)

        !if denominator not zero
        IF ( ABS(dp_dwn).GT.1.0E-20 ) THEN
            lmbd = dp_up / dp_dwn
        END IF

        calcLambda = lmbd
      
      END FUNCTION calcLambda

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE PartInElement2(x,centers,normals,nf,res,ierr)
!
!      
! -----------------------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER ierr,nf
      REAL(8) x(3), centers(nf,3), normals(nf,3), res, maxres

      INTEGER i

      ierr = 0
      maxres = -1

      DO i=1,nf

        CALL DotProduct(normals(i,:),x-centers(i,:),res)
        !WRITE(*,'(A5,F14.6)') "res=", res
        IF (res.GT.0.0D0) ierr = 1

        IF (res.GT.maxres) maxres = res

      END DO

      res = maxres

  

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE CalcTrackLambda(a,b,Cf,norm,lambda,ierr)
!
!     Calculates the ratio of tracking step ("track") from start to the 
!     intersection with a face (its plane actually)
!     a = starting point vector (x,y,z)
!     b = end point vector (x,y,z)
!     Cf = face center vector (x,y,z)
!     norm = face normal vector (x,y,z)
!     lambda = resulting ratio; valid = 0..1; invalid (no intersection) = <0 or >1
!     ierr = error flag 
!            
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE 
      INTEGER ierr
      REAL(8) a(3),b(3),Cf(3),norm(3),lambda

      REAL(8) dp_up, dp_dwn

      lambda = -1.0D0

      

      !numerator
      CALL DotProduct(Cf-a,norm,dp_up)

      !denominator
      CALL DotProduct(b-a,norm,dp_dwn)
      

!      PRINT *, "calc tracking lambda prints:"
!      WRITE(*,'(A6,3G12.4)') " a=",a
!      WRITE(*,'(A6,3G12.4)') " b=",b
!      WRITE(*,'(A6,3G12.4)') "b-a=",b-a
!      WRITE(*,'(A6,3G12.4)') " Cf=",Cf
!      WRITE(*,'(A6,3G12.4)') "Cf-a=",Cf-a
!      WRITE(*,'(A6,3G12.4)') " n=",norm
!      WRITE(*,'(A6,G12.4)') " dpU=",dp_up
!      WRITE(*,'(A6,G12.4)') " dpD=",dp_dwn
!      PRINT *, "----------------------------"



      !if denominator not zero
      IF ( ABS(dp_dwn).GT.1.0E-20 ) THEN

            lambda = dp_up / dp_dwn
            ierr = 0

            RETURN
      ELSE
            PRINT *, "ERROR :: CalcTrackLambda :: Denominator is zero!"
      END IF

      ierr = 1


      END