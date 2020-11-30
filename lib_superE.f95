! -----------------------------------------------------------------------------------------
      SUBROUTINE GetFlowField(m,fluid)
!
!     Get fluid flow field
!
! -----------------------------------------------------------------------------------------
      USE inpFile
      USE mesh
      USE mFluid   

      IMPLICIT NONE

      TYPE(meshType) :: m  ! mesh data structure
      TYPE(fluidType) :: fluid ! fluid results    
      
      INTEGER i

      IF (inp%AnalFlowField.EQ.0) THEN 
!
!       Read Fluid flow field
!
        CALL ReadVTKfluid(m,fluid)
!
!       Calculate flow gradients
!          
        CALL SetGrad2(m,fluid%Un(:,1),fluid%gradUxn)
        CALL SetGrad2(m,fluid%Un(:,2),fluid%gradUyn)
        CALL SetGrad2(m,fluid%Un(:,3),fluid%gradUzn)  

        !CALL gradTesting(m,fluid)

      ELSE
!
!       Produce analytical flow field
!      
        CALL ProduceAnalFF(m,fluid)
      END IF
!
!     Get vorticity from velocity gradient
!      
      IF (fluid%iVortn.EQ.0) THEN
        CALL GetWfromVgrad(fluid)
        fluid%iVortn=1
      END IF

!      norm = 0.0D0
!      norm2 = 0.0D0
!      Do i=1,VTKmesh%nnodes
!        norm = norm + ( fluid%gradUxn(i,1) - 2.0*(0.5D0+VTKmesh%x(i,1))   )**2.0D0
!        norm2 = norm2 + ( fluid%gradUxn(i,1) ) ** 2.0D0
!      end do
!      print *,SQRT(norm/norm2)

!
!     Get polyhedron volumetric center values
!
      !CALL PolyVolCenterVals(m,fluid)

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE gradTesting(m,fluid)
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE mFluid
      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure
      TYPE(fluidType) :: fluid ! fluid results

      INTEGER i,kx,ky,kz

      REAL(8) x,y,z,w,pi,A,B,C
      REAL(8) dudx,dudy,dudz
      REAL(8) dvdx,dvdy,dvdz
      REAL(8) errX,errY,errZ

      w = 1.0D0
      pi = 3.14159D0
      A = -1.02D0
      B = 1.01D0
      C = 0.98D0

      !CALL SetGrad(m,fluid%Un(:,1),fluid%gradUxn)

      errX = 0.0D0
      errY = 0.0D0
      errZ = 0.0D0

      kx=0
      ky=0
      kz=0
      DO i = 1,m%nnodes

        x = m%x(i,1)
        y = m%x(i,2)
        z = m%x(i,3)
        
        dudx = A*w*pi*cos(w*pi*(A*x+B*y+C*z))
        dudy = B*w*pi*cos(w*pi*(A*x+B*y+C*z))
        dudz = C*w*pi*cos(w*pi*(A*x+B*y+C*z))

        dvdx = fluid%gradUxn(i,1)
        dvdy = fluid%gradUxn(i,2)
        dvdz = fluid%gradUxn(i,3)

        IF (dudx .NE. 0) THEN
          kx = kx + 1
          errX = errX + ABS( (dvdx - dudx) / dudx )
        END IF

        IF (dudy .NE. 0) THEN
          ky = ky + 1
          errY = errY + ABS( (dvdy - dudy) / dudy )
        END IF

        IF (dudz .NE. 0) THEN
          kz = kz + 1
          errZ = errZ + ABS( (dvdz - dudz) / dudz )
        END IF

        !WRITE(*,'(A5,F16.8,A5,F12.5,A5,F12.5)') " ana=",dudx," dis=",dvdx

      END DO

      errX = errX / kx 
      errY = errY / ky 
      errZ = errZ / kz 

      WRITE(*,'(A10,F12.6)') "avgerrX=",errX
      WRITE(*,'(A10,F12.6)') "avgerrY=",errY
      WRITE(*,'(A10,F12.6)') "avgerrZ=",errZ

      PRINT *, "Stopping after grad test"
      CALL stopProgram(1)


      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE AllocateFluidFileds(m,fluid)
!
!     Allocate memory for flow fields
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE mFluid
      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure
      TYPE(fluidType) :: fluid ! fluid results

      fluid%nelem=m%nelem
      fluid%nnodes=m%nnodes
      ALLOCATE (fluid%Un(fluid%nnodes,3))
      ALLOCATE (fluid%Vortn(fluid%nnodes,3))
      ALLOCATE (fluid%Ue(fluid%nelem,3))
      ALLOCATE (fluid%Pn(fluid%nnodes))
      ALLOCATE (fluid%Tn(fluid%nnodes))
      ALLOCATE (fluid%Pe(fluid%nelem))
      ALLOCATE (fluid%gradUxn(fluid%nnodes,3))
      ALLOCATE (fluid%gradUyn(fluid%nnodes,3))
      ALLOCATE (fluid%gradUzn(fluid%nnodes,3)) 
      ALLOCATE (fluid%dvxdt(fluid%nnodes))
      ALLOCATE (fluid%dvydt(fluid%nnodes))
      ALLOCATE (fluid%dvzdt(fluid%nnodes))
!
!     Init fluid fields
!
      fluid%Un=0.0D0
      fluid%Ue=0.0D0
      fluid%Pn=0.0D0            
      fluid%Pe=0.0D0     
      fluid%gradUxn=0.0D0 
      fluid%gradUyn=0.0D0
      fluid%gradUzn=0.0D0
      fluid%Vortn=0.0D0
      fluid%dvxdt=0.0D0
      fluid%dvydt=0.0D0
      fluid%dvzdt=0.0D0
!
!     Set switches
!
      fluid%iUn=0
      fluid%iUe=0
      fluid%iPn=0
      fluid%iPe=0
      fluid%iTn=0
      fluid%iVortn=0

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE GetWfromVgrad(fluid)
!
!     Get vorticity from velocity gradient
!      
! -----------------------------------------------------------------------------------------
      USE mFluid
      IMPLICIT NONE 

      INTEGER i
      TYPE(fluidType) :: fluid ! fluid results

      DO i=1,fluid%nnodes
        fluid%Vortn(i,1) = fluid%gradUzn(i,2)  - fluid%gradUyn(i,3)
        fluid%Vortn(i,2) = fluid%gradUxn(i,3)  - fluid%gradUzn(i,1)
        fluid%Vortn(i,3) = fluid%gradUyn(i,1)  - fluid%gradUxn(i,2)
      END DO

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE ProduceAnalFF(m,fluid)
!
!     Allocate memory for flow fields
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE mFluid
      USE inpFile      
      USE logFile            
      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure
      TYPE(fluidType) :: fluid ! fluid results
      REAL(8) x,y,z,r
      INTEGER i

      CALL logWrite ("MESSAGE :: ProduceAnalFF :: setting up analytical flow field!")
 
      IF (inp%AnalFlowField.EQ.1) THEN 

        CALL logWrite ("MESSAGE :: ProduceAnalFF :: u=(x,y,-2z)")
        fluid%iUn=1
        DO i=1,m%nnodes
          x=m%x(i,1)
          y=m%x(i,2)
          z=m%x(i,3)
          fluid%Un(i,1)=x
          fluid%Un(i,2)=y
          fluid%Un(i,3)=-2.0D0*z
          fluid%gradUxn(i,1)=1.0D0 
          fluid%gradUxn(i,2)=0.0D0 
          fluid%gradUxn(i,3)=0.0D0 

          fluid%gradUyn(i,1)=0.0D0
          fluid%gradUyn(i,2)=1.0D0
          fluid%gradUyn(i,3)=0.0D0
    
          fluid%gradUzn(i,1)=0.0D0
          fluid%gradUzn(i,2)=0.0D0
          fluid%gradUzn(i,3)=-2.0D0
        END DO

      ELSE IF (inp%AnalFlowField.EQ.2) THEN ! Couette

        CALL logWrite ("MESSAGE :: ProduceAnalFF :: u=(z,0,0)")
        fluid%iUn=1
        DO i=1,m%nnodes
          x=m%x(i,1)
          y=m%x(i,2)
          z=m%x(i,3)
          fluid%Un(i,1)=z
          fluid%Un(i,2)=0.0D0
          fluid%Un(i,3)=0.0D0
          fluid%gradUxn(i,1)=0.0D0 
          fluid%gradUxn(i,2)=0.0D0 
          fluid%gradUxn(i,3)=1.0D0 

          fluid%gradUyn(i,1)=0.0D0
          fluid%gradUyn(i,2)=0.0D0
          fluid%gradUyn(i,3)=0.0D0
    
          fluid%gradUzn(i,1)=0.0D0
          fluid%gradUzn(i,2)=0.0D0
          fluid%gradUzn(i,3)=0.0D0
        END DO

      ELSE IF (inp%AnalFlowField.EQ.3) THEN ! Pipe flow

        CALL logWrite ("MESSAGE :: ProduceAnalFF :: u=(0,0,8*(1/4-r^2))")
        fluid%iUn=1
        DO i=1,m%nnodes
          x=m%x(i,1)
          y=m%x(i,2)
          z=m%x(i,3)
          fluid%Un(i,1)=0.0D0
          fluid%Un(i,2)=0.0D0
          fluid%Un(i,3)=8.0D0*(0.25D0-(x*x+y*y))
          fluid%gradUxn(i,1)=0.0D0 
          fluid%gradUxn(i,2)=0.0D0 
          fluid%gradUxn(i,3)=0.0D0 

          fluid%gradUyn(i,1)=0.0D0
          fluid%gradUyn(i,2)=0.0D0
          fluid%gradUyn(i,3)=0.0D0
    
          fluid%gradUzn(i,1)=-16.0D0*x
          fluid%gradUzn(i,2)=-16.0D0*y
          fluid%gradUzn(i,3)=0.0D0
        END DO

        ELSE IF (inp%AnalFlowField.EQ.4) THEN ! uniform in z direction

        CALL logWrite ("MESSAGE :: ProduceAnalFF :: u=(0,0,1)")
        fluid%iUn=1
        DO i=1,m%nnodes
          x=m%x(i,1)
          y=m%x(i,2)
          z=m%x(i,3)
          fluid%Un(i,1)=0.0D0
          fluid%Un(i,2)=0.0D0
          fluid%Un(i,3)=1.0D0
          fluid%gradUxn(i,1)=0.0D0 
          fluid%gradUxn(i,2)=0.0D0 
          fluid%gradUxn(i,3)=0.0D0 

          fluid%gradUyn(i,1)=0.0D0
          fluid%gradUyn(i,2)=0.0D0
          fluid%gradUyn(i,3)=0.0D0
    
          fluid%gradUzn(i,1)=0.0D0
          fluid%gradUzn(i,2)=0.0D0
          fluid%gradUzn(i,3)=0.0D0
        END DO

        ELSE IF (inp%AnalFlowField.EQ.5) THEN ! Couette

        CALL logWrite ("MESSAGE :: ProduceAnalFF :: u=(0,0,x)")
        fluid%iUn=1
        DO i=1,m%nnodes
          x=m%x(i,1)
          y=m%x(i,2)
          z=m%x(i,3)
          fluid%Un(i,1)=0.0D0
          fluid%Un(i,2)=0.0D0
          fluid%Un(i,3)=x+0.5D0
          fluid%gradUxn(i,1)=0.0D0 
          fluid%gradUxn(i,2)=0.0D0 
          fluid%gradUxn(i,3)=0.0D0 

          fluid%gradUyn(i,1)=0.0D0
          fluid%gradUyn(i,2)=0.0D0
          fluid%gradUyn(i,3)=0.0D0
    
          fluid%gradUzn(i,1)=1.0D0
          fluid%gradUzn(i,2)=0.0D0
          fluid%gradUzn(i,3)=0.0D0
        END DO

      ELSE
        CALL logWrite ("ERROR :: ProduceAnalFF :: unknown AnalFlowField")
        CALL StopProgram(1)
      END IF

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE PeriodicBC(part)
!
!     Write to BIN file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE counters 
      USE logFile   
      USE ParalelEnvironment     
      IMPLICIT NONE

      TYPE(SuperElType) part

      IF (inp%Periodic.EQ.1) THEN
        IF (part%r(inp%PeriodicAxis).GT.inp%PeriodicFrom) THEN
          part%r(inp%PeriodicAxis) = part%r(inp%PeriodicAxis) -  ( inp%PeriodicFrom - inp%PeriodicTo )
        END IF
      ELSE IF (inp%Periodic.EQ.2) THEN
        IF (part%r(inp%PeriodicAxis).LT.inp%PeriodicFrom) THEN
          part%r(inp%PeriodicAxis) = part%r(inp%PeriodicAxis) -  ( inp%PeriodicFrom - inp%PeriodicTo )
        END IF
      ELSE
        CALL logIntWrite("ERROR :: PeriodicBC :: Wrong code : ",inp%Periodic)
        CALL StopProgram(1)
      END IF

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE ExportToBIN(part)
!
!     Write to BIN file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE counters 
      USE logFile   
      USE ParalelEnvironment     
      IMPLICIT NONE

      TYPE(SuperElType) part(inp%maxNp)
      INTEGER lun,i
      CHARACTER(255) vrstica

      lun = 12      
      WRITE(vrstica,'(A,A,A,A,I0,A,I0,A)') TRIM(inp%BINResultsFolder),"/",TRIM(lg%IDname),"-",cnt%iTime,"-",env%mpr,".bin"   
      OPEN (lun,FILE=TRIM(vrstica),ERR=10,FORM='UNFORMATTED',STATUS='UNKNOWN')

      WRITE (lun) cnt%cnp
      DO i=1,cnt%cnp
            CALL seWriteBIN(lun,part(i))
      END DO

      CLOSE (lun)

      RETURN

10    CONTINUE ! error when opening input file
      CALL logWrite ("ERROR :: ExportToBIN :: Could not open file for export!")
      CALL StopProgram(1)

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE ExportResultsToASCII(part,m,fluid)
!
!     Write to file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mesh
      USE logFile      
      USE mFluid   
      USE counters   
      USE ParalelEnvironment      
      IMPLICIT NONE

      TYPE(meshType) :: m
      TYPE(fluidType) :: fluid ! fluid results
      TYPE(SuperElType) part(inp%maxNp)

      INTEGER iT,lun,mpr,cnp,i
      CHARACTER(255) vrstica,vrstica2
      LOGICAL exists

      CALL logWrite("Exporting results to ASCII files!")
!
!     Loop over BIN files and export
!
      lun = 12  

      cnt%iTime=-1
      DO iT=0,inp%nTimeSteps
        cnt%iTime = cnt%iTime + 1
        IF (cntExport(inp%ResultExportFreq)) THEN
!
!         Check if all proc files exist
!
            DO mpr=0,env%nproc-1
                  WRITE(vrstica,'(A,A,A,A,I0,A,I0,A)') TRIM(inp%BINResultsFolder),"/",TRIM(lg%IDname),"-",cnt%iTime,"-",mpr,".bin"
                  INQUIRE(FILE=TRIM(vrstica),EXIST=exists)
                  IF(.NOT.exists) THEN
                        !PRINT *, "ExportResultsToASCII :: Skipping timestep",iT
                        GOTO 9
                  END IF
            END DO
!
!         Read a single time step
!
          cnt%cnp=0
          DO mpr=0,env%nproc-1
            WRITE(vrstica,'(A,A,A,A,I0,A,I0,A)') TRIM(inp%BINResultsFolder),"/",TRIM(lg%IDname),"-",cnt%iTime,"-",mpr,".bin"   
            OPEN (lun,FILE=TRIM(vrstica),ERR=10,FORM='UNFORMATTED',STATUS='OLD')

                  READ (lun) cnp
                  DO i=1,cnp
                        CALL seReadBIN(lun,part(cnt%cnp+i))  ! read particle from binary file
                        CALL seBaseInit(part(cnt%cnp+i))     ! calculate info not in binary file
                  END DO
                  cnt%cnp = cnt%cnp + cnp            

                  CLOSE (lun)  

          END DO
!
!         Export time step to ascii files
!
          CALL ExportResults(part,m,fluid)
        END IF

9     CONTINUE

      END DO

      RETURN

10    CONTINUE ! error when opening input file
      WRITE(vrstica2,'(A,A)') "ERROR :: ExportResultsToASCII :: Could not open BIN file: ",TRIM(vrstica)
!      !CALL logWrite ("ERROR :: ExportResultsToASCII :: Could not open BIN file!")
      CALL logWrite (TRIM(vrstica2))
!      PRINT *, "ls"
!      CALL SYSTEM("ls")
!      PRINT *, "ls binary"
!      CALL SYSTEM("ls binary")
      CALL StopProgram(1)
!

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE ExportResults(part,m,fluid)
!
!     Write to file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mesh
      USE logFile      
      USE mFluid   
      USE counters         
      IMPLICIT NONE

      TYPE(meshType) :: m
      TYPE(fluidType) :: fluid ! fluid results
      TYPE(SuperElType) part(inp%maxNp)

      IF (inp%TecplotMeshExport.GT.0) THEN
        IF (MODULO(cnt%iTime,inp%TecplotMeshExport).EQ.0) CALL TecplotMeshExport(part)
      END IF

      IF (inp%ParaviewMeshExport.GT.0) THEN
        IF (MODULO(cnt%iTime,inp%ParaviewMeshExport).EQ.0) CALL ParaviewMeshExport(part)
      END IF   
  
      IF (inp%ParaviewListExport.GT.0) THEN
        IF (MODULO(cnt%iTime,inp%ParaviewListExport).EQ.0) CALL ParticleListExportCSV(part)
      END IF

      IF (inp%ParaviewListVTKExport.GT.0) THEN
        IF (MODULO(cnt%iTime,inp%ParaviewListVTKExport).EQ.0) CALL ParticleListExportVTK(part)
      END IF

      IF (inp%WriteVTKflowResults.GT.0) THEN
        IF (MODULO(cnt%iTime,inp%WriteVTKflowResults).EQ.0) CALL WriteVTKflowResults(m,fluid)
      END IF

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE TecplotMeshExport(part)
!
!     Write to file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE logFile
      USE mesh
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) part(inp%maxNp)
      TYPE(localMeshType)  :: sphere
      REAL(8) time,r(3)
      INTEGER np,i,j,k,lun
      CHARACTER(255) vrstica


      CALL MeshAsphere(sphere)
  

      lun = 22
!
!     Count number of particles for export
!
      np=cnt%cnp
!
!     First time header
!
      WRITE(vrstica,'(A,A)') TRIM(lg%IDname),".tec-mesh.dat"
      IF (inp%iTME.EQ.1) THEN    
            OPEN (lun,FILE=TRIM(vrstica),ERR=10,STATUS='UNKNOWN')
            WRITE (lun,'(A)') 'VARIABLES = "X", "Y", "Z"'
!           write time step header
            WRITE (lun,'(A,A,F10.5,A,I8,A,I8)') &
     &        'ZONE T="SuperE Particles",', &
     &        'StrandID=1, SOLUTIONTIME=',cnt%rTime, &
     &        ', F=FEPOINT, ET=TRIANGLE,N= ', &
             np*sphere%nnodes, ',E=', np*sphere%nbelem
      ELSE
        OPEN (lun,FILE=TRIM(vrstica),ERR=10,ACCESS='APPEND',STATUS='OLD')
!       write time step header
        WRITE (lun,'(A,A,F10.5,A,I8,A,I8,A)') &
     &        'ZONE T="Elliptic Particles",', &
     &        'StrandID=1, SOLUTIONTIME=',cnt%rTime, &
     &        ', F=FEPOINT, ET=TRIANGLE,N= ', &
             np*sphere%nnodes, ',E=', np*sphere%nbelem, &
!     &      ', VARSHARELIST=([3]=1), CONNECTIVITYSHAREZONE = 1' &
           ', CONNECTIVITYSHAREZONE = 1'
      END IF


!     Loop over mesh nodes
      DO i=1,np
          DO j=1,sphere%nnodes
!           Transform shpere to ellipsoid (diameter of sphere = 1)
            CALL seMapSphereToSuperE(part(i),r,2.0D0*sphere%x(j,:))
!           Rotate and move to particle location
            CALL seTransfL2G(part(i),r)

            WRITE (vrstica,*) r(1),r(2),r(3)
            CALL sqblnk(lun,vrstica)
          END DO
      END DO


      IF (inp%iTME.EQ.1) THEN 
        DO i=1,np
            DO j=1,sphere%nbelem
              WRITE (vrstica,*) ( (i-1)*sphere%nnodes+sphere%ibc(j,k),k=1,sphere%npob)
              CALL sqblnk(lun,vrstica)
            END DO
        END DO
        inp%iTME=0        
      END IF

      CLOSE (lun)

      RETURN

10    CONTINUE ! error when opening input file
      CALL logWrite ("ERROR :: TecplotMeshExport :: Could not open file for export!")
      CALL StopProgram(1)

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE ParaviewMeshExport(part)
!
!     Write to file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE logFile
      USE mesh
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) part(inp%maxNp)
      TYPE(localMeshType)  :: sphere
      REAL(8) time,r(3)
      INTEGER np,i,j,k,lun
      CHARACTER(255) vrstica
      INTEGER iTime

      iTime=1
!
!     Create subdirectory for Paraview files
!
!      CALL create_directory( inp%ParaviewResultsFolder )

      CALL MeshAsphere(sphere)
  

      lun = 22
!
!     Count number of particles for export
!
      np=cnt%cnp
!
!     First time header
!
      WRITE(vrstica,'(A,A,A,A,I0,A)') TRIM(inp%ParaviewResultsFolder),"/",TRIM(lg%IDname),".para-mesh-",cnt%iTime,".vtk"   
      OPEN (lun,FILE=TRIM(vrstica),ERR=10,STATUS='UNKNOWN')
      WRITE (lun,'(A)') '# vtk DataFile Version 2.0'
      WRITE (lun,'(A)') 'Superellipsoid mesh plot'
      WRITE (lun,'(A)') 'ASCII'
      WRITE (lun,'(A)') 'DATASET UNSTRUCTURED_GRID'
      WRITE (lun,'(A,I10,A)') 'POINTS',np*sphere%nnodes," float"

!     Loop over mesh nodes
      DO i=1,np
          DO j=1,sphere%nnodes
!           Transform shpere to ellipsoid (diameter of sphere = 1)
            CALL seMapSphereToSuperE(part(i),r,2.0D0*sphere%x(j,:))
!           Rotate and move to particle location
            CALL seTransfL2G(part(i),r)

            WRITE (vrstica,*) r(1),r(2),r(3)
            CALL sqblnk(lun,vrstica)
          END DO
      END DO

      WRITE (lun,'(A,I10,I10)') 'CELLS ',np*sphere%nbelem,np*sphere%nbelem*4
 
        DO i=1,np
            DO j=1,sphere%nbelem
              WRITE (vrstica,*) "3 ",( (i-1)*sphere%nnodes+sphere%ibc(j,k)-1,k=1,sphere%npob)
              CALL sqblnk(lun,vrstica)
            END DO
        END DO

      WRITE (lun,'(A,I10)') 'CELL_TYPES ',np*sphere%nbelem

        DO i=1,np
            DO j=1,sphere%nbelem
              WRITE (vrstica,*) "5"
              CALL sqblnk(lun,vrstica)
            END DO
        END DO

      CLOSE (lun)

      RETURN

10    CONTINUE ! error when opening input file
      CALL logWrite ("ERROR :: ParaviewMeshExport :: Could not open file for export!")
      CALL StopProgram(1)

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE seWriteToFile(se)
!
!     Write to file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE mesh
      IMPLICIT NONE

      TYPE(SuperElType) se
      TYPE(localMeshType)  :: sphere

      REAL(8), ALLOCATABLE :: RM(:,:), RMTrans(:,:) ! rotation matrices
      REAL(8), ALLOCATABLE :: r(:) ! a vector
      REAL(8) time
      CHARACTER(255) vrstica
      INTEGER i,j,k,lun,ip,np

      ALLOCATE (RM(3,3),RMTrans(3,3),r(3))

      time=3.3D0
      np=1
      lun=99

!     create spherical mesh ( D=1, r=0.5, c=(0,0,0) )
      CALL MeshAsphere(sphere)
!
!     Particle header
!
!     Write file header
      OPEN (lun,FILE=TRIM("tri.se.dat"),STATUS='UNKNOWN')
      WRITE (lun,'(A)') 'VARIABLES = "X", "Y", "Z"'
!     write time step header
      WRITE (lun,'(A,A,F10.5,A,I8,A,I8)') &
     &        'ZONE T="SuperE Particles",', &
     &        'StrandID=1, SOLUTIONTIME=',time, &
     &        ', F=FEPOINT, ET=TRIANGLE,N= ', &
             np*sphere%nnodes, ',E=', np*sphere%nbelem

!
!     Export node list
!
!     Calculate rotation matrix
!      CALL seCalRotationMatrix(se) ! has alerady been calculated

!     Loop over mesh nodes
          DO j=1,sphere%nnodes
!           Transform shpere to ellipsoid (diameter of sphere = 1)
            CALL seMapSphereToSuperE(se,r,2.0D0*sphere%x(j,:))
!           Rotate and move to particle location
            CALL seTransfL2G(se,r)

            WRITE (vrstica,*) r(1),r(2),r(3)
            CALL sqblnk(lun,vrstica)
          END DO

!
!     Export particle mesh connectivity (first time only)
!
!       For each particle
        ip=0
            DO j=1,sphere%nbelem
              WRITE (vrstica,*) (ip*sphere%nnodes+sphere%ibc(j,k),k=1,sphere%npob)
              CALL sqblnk(lun,vrstica)
            END DO


!     Close file
      CLOSE(lun)


      DEALLOCATE (RM,RMTrans,r)

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE AddParticles(m,part,fluid)
!
!     Add particles
!      
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mesh
      USE logFile
      USE mFluid
      USE counters
      USE ParalelEnvironment
      IMPLICIT NONE

      TYPE(meshType) :: m  ! mesh data structure
      TYPE(SuperElType) part(inp%maxNp)
      TYPE(fluidType) :: fluid
      INTEGER n,model,when,i,j,ierr
      REAL(8) a,b,c,e1,e2,rho,box(3,2),u(3)
      LOGICAL CheckKdaj
      CHARACTER*255 patchName
      CHARACTER*255 uselesstext

      j=0
      ierr = 1

      DO i=1,inp%nPartStat
        READ(inp%PartStat(i),*) model,when
        IF ( CheckKdaj(cnt%iTime,when).EQV..TRUE.) THEN
!
!         Add a particle
!

!         Model = 1 -> AddParticlesFromPatch
!         Add particles positined randomly accros the patch
          IF (model.EQ.1) THEN

            READ(inp%PartStat(i),*) model,when,n,a,b,c,e1,e2,rho,u,patchName

            CALL AddParticlesFromPatch(m,fluid,part,n,a,b,c,e1,e2,rho,u,patchName,ierr) 
            IF (ierr.EQ.0) j = j+n
            
!         Model = 2 -> AddParticlesFromPatch -limit location by geometric shape
!         Add particles positined randomly accros the limited space
          ELSE IF (model.EQ.2) THEN

            !READ(inp%PartStat(i),*) model,when,n,a,b,c,e1,e2,rho,u,patchName
            !CALL AddParticlesFromPatch2(m,fluid,part,n,a,b,c,e1,e2,rho,u,patchName,ierr) 

            READ(inp%PartStat(i),*) model,when,n,a,b,c,e1,e2,rho,u,patchName
            CALL AddParticlesFromPatch2(m,fluid,part,n,a,b,c,e1,e2,rho,u,patchName,ierr)
            IF (ierr.EQ.0) j = j+n


!         Model = 3 -> AddParticleKnownPosition
!         Set one by one particle, by finding nearest element to particle location
!         and searching first in this element, and following neighbour elements
          ELSE IF (model.EQ.3) THEN

            READ(inp%PartStat(i),*) model,when,n,a,b,c,e1,e2,rho &
     &                             ,box(1,1),box(1,2),box(2,1)

            CALL AddParticleKnownPosition(m,fluid,part,n,a,b,c,e1,e2,rho,box,ierr) 
            IF (ierr.EQ.0) j = j+n

          ELSE
            CALL logIntWrite ("ERROR :: AddParticles :: Unknown model : ",model)
            CALL StopProgram(1)
          END IF
        END IF
      END DO

      IF (j.GT.0) WRITE(*,'(A23,I6)') "Tot. particles added::",j

      !some init stuff
      DO i = 1,j
        part(i)%nnf = 0
      END DO
      

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE EstimateNumberOfParticles()
!
!     Based on input file, estimate the total number of particles
!      
! -----------------------------------------------------------------------------------------
      USE inpFile
      USE logFile
      USE counters
      IMPLICIT NONE

      INTEGER n,model,when,i,iT
      LOGICAL CheckKdaj

      inp%maxNp = 0

      DO iT=0,inp%nTimeSteps
        DO i=1,inp%nPartStat
          READ(inp%PartStat(i),*) model,when
          IF ( CheckKdaj(iT,when).EQV..TRUE.) THEN

            IF (model.GT.0.AND.model.LE.3) THEN

              READ(inp%PartStat(i),*) model,when,n
              inp%maxNp = inp%maxNp + n
            ELSE
              CALL logIntWrite ("ERROR :: EstimateNumberOfParticles :: Unknown model : ",model)
              CALL StopProgram(1)
            END IF
          END IF
        END DO
      END DO

      END




! -----------------------------------------------------------------------------
      LOGICAL FUNCTION CheckKdaj(tstep,kdaj)
      INTEGER kdaj,tstep
      LOGICAL add
      
      add=.FALSE.
      IF (kdaj.EQ.0.AND.tstep.EQ.0) THEN
        add=.TRUE.
      END IF
      IF (kdaj.GT.0) THEN
        IF (MOD(tstep,kdaj).EQ.0) THEN
         add=.TRUE.
       END IF
      END IF
      
      CheckKdaj=add
      
      RETURN
      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE AddParticlesFromPatch2(m,fluid,part,n,a,b,c,e1,e2,rho,u,patchName,ierr)
!
!     Add random particles
!      
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mesh
      USE logFile
      USE counters
      USE ParalelEnvironment
      USE mFluid

      IMPLICIT NONE
      
      TYPE(meshType) :: m  ! mesh data structure
      TYPE(SuperElType) part(inp%maxNp),p
      TYPE(fluidType) :: fluid
      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position   

      INTEGER i,j,k,n,ierr
      REAL(8) a,b,c,e1,e2,rho,u(3)
      CHARACTER*64 patchName


      INTEGER ibound,sf,ef
      INTEGER ie,iface, itri

      REAL(8) dist,xp(3),xb(3),xtri(3,3),xcorr(3)

      REAL(8) xr,xalpha

      REAL(8) mins(3),maxs(3),dx,avg_vv(3),len

      
      REAL(8) r(3),v(3),o(3),x,phi,theta,psi,pi

      DO i = 1,m%nbound 

        !PRINT *, "bound(",i,")=",m%bound(i)%name 
        IF ( TRIM(m%bound(i)%name).EQ.TRIM(patchName) ) ibound = i

      END DO

      PRINT *, "Found patch at i =",i

      !Start and end face
      sf = m%bound(ibound)%startFace + 1
      ef = sf + m%bound(ibound)%nFaces - 1


      !Get patch bounding box
      k=0
      mins = 9E+9
      maxs = -9E+9
      avg_vv = 0.0D0
      DO i = sf,ef

        DO j = 1,3
          IF(m%f(i)%center(j).LT.mins(j)) mins(j) = m%f(i)%center(j)
          IF(m%f(i)%center(j).GT.maxs(j)) maxs(j) = m%f(i)%center(j)
        END DO

        avg_vv = avg_vv + fluid%Un( m%f(i)%con(1),: )
        k = k + 1

      END DO

      avg_vv = avg_vv / k
      CALL vecLen(avg_vv,len)

      avg_vv = avg_vv / len

      mins = mins + 1E-06*avg_vv
      maxs = maxs + 1E-06*avg_vv


      PRINT *, "Domain span"
      WRITE(*, '(A4,3F12.6)') "min=",mins
      WRITE(*, '(A4,3F12.6)') "max=",maxs


      !Loop for number of particles 
      k = 0
      DO i = 1,n

        p = part(i)

10      CONTINUE 
        k = k + 1

        IF(k.GT.5*n) GOTO 11

        CALL RANDOM_NUMBER(x)
        dx = maxs(1) - mins(1)
        p%r(1) = mins(1) + x*dx

        CALL RANDOM_NUMBER(x)
        dx = maxs(2) - mins(2)
        p%r(2) = mins(2) + x*dx

        CALL RANDOM_NUMBER(x)
        dx = maxs(3) - mins(3)
        p%r(3) = mins(3) + x*dx

        !Find particle across boundary faces
        DO j = sf,ef 

          ie = FindInFace(p%r,j,xcorr)   

          !If found exit loop
          IF (ie.GT.0) THEN 
            p%r = xcorr
            EXIT
          END IF

        END DO

        IF (ie.LT.1) GOTO 10

        p%element = ie
        !PRINT *, "Particle r:", p%r
        !PRINT *,  "Found at:",ie

        !Try if found, else go back, and try again with another one
        CALL GetFFFap(m,p,fluid,fap,-1,ierr)
        IF (ierr.NE.0) THEN 
          !PRINT *, "ierr ne 0, trying again"
          GOTO 10
        ELSE
          !WRITE(*,'(A,3F10.5)') "Found at interpol, fapvx=",fap%vx
        END IF

        !Set particle orientation
        v=0.0D0
        pi=4.0D0*ATAN(1.0D0)
        CALL RANDOM_NUMBER(x)
        theta=ACOS(2.0D0*x-1.0D0)  ! -pi .. pi
        CALL RANDOM_NUMBER(x)
        phi=2.0D0*pi*x ! ni fajn ce je nic
        CALL RANDOM_NUMBER(x)
        psi=2.0D0*pi*x
        o=0.0D0

        r = p%r
        CALL seInit2(p,a,b,c,e1,e2,rho,phi,theta,psi,r,v,o)
        !PRINT *, "after seinit"
        !PRINT *, "p%e=",p%element
     
        !If parsed nonzero, set parsed velocity
        IF ( u(1).NE.0.0D0.OR.u(2).NE.0.0D0.OR.u(3).NE.0.0D0 ) THEN
          p%v(1)=u(1)
          p%v(2)=u(2)
          p%v(3)=u(3)
        
        !If parsed zero, set to flow velocity
        ELSE
          CALL GetFFFap(m,p,fluid,fap,-1,ierr)
          p%v(1)=fap%vx
          p%v(2)=fap%vy
          p%v(3)=fap%vz
    
          p%vOld = p%v
    
          !WRITE(*,'(A4,3F10.5)') "p%r=",p%r
          !WRITE(*,'(A4,3F10.5)') "p%v=",p%v
        END IF
  
        p%pid=par_NextPID(cnt%LastPID)
        p%id=p%pid*100000000+cnt%cnp
        cnt%LastPID = p%pid
        
        IF (par_IsMine(p%pid)) THEN 
          cnt%cnp=cnt%cnp+1
          part(cnt%cnp)=p
        END IF

      END DO

      ierr = 0
      RETURN

11    CONTINUE

      ierr = 0
      PRINT *, "Failed to insert all particles"
      PRINT *, "N parts =", i, " percentage=", INT( 100*i/n )
      PRINT *, "It's possible to continue if you want, with less parts"

      CALL BreakPoint()

      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      INTEGER FUNCTION FindInFace(x,iface,xnew)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: iface
        REAL(8), INTENT(IN) :: x(3)
        REAL(8), INTENT(OUT) :: xnew(3)

        INTEGER i,ie,ntri,ierr

        REAL(8) xtri(3,3),fc(3),fn(3),fb(3)

        REAL(8) dp,frac 

        frac = 2.5E-02

        FindInFace = -1
        ntri = size( m%f(iface)%tf )

        DO i = 1,ntri

          xtri(1,:) = m%x( m%f(iface)%tf(i)%con(1),: )
          xtri(2,:) = m%x( m%f(iface)%tf(i)%con(2),: )
          xtri(3,:) = m%x( m%f(iface)%tf(i)%con(3),: )

          fc = m%f(iface)%tf(i)%center
          fn = m%f(iface)%tf(i)%normal
    
          CALL GetBarycentricInTri(x,xtri,fc,fn,fb,ierr)

          IF(ierr.EQ.0) THEN

            ie = m%f(iface)%owner

            xnew = fb(1)*xtri(1,:) + fb(2)*xtri(2,:) + fb(3)*xtri(3,:)

            !Move particle by 'frac' of the perpendicular distance from 
            !location on the face towards the element center, to prevent
            !boundary hits at injection
            CALL DotProduct(fn, m%e(ie)%xc-fc, dp)
            xnew = xnew + frac*fn*dp

            FindInFace = ie
            RETURN

          END IF

        END DO
    
      END FUNCTION FindInFace

      END 



! -----------------------------------------------------------------------------------------
      SUBROUTINE AddParticlesFromPatch(m,fluid,part,n,a,b,c,e1,e2,rho,u,patchName,ierr)
!
!     Add random particles
!      
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mesh
      USE logFile
      USE counters
      USE ParalelEnvironment
      USE mFluid

      IMPLICIT NONE
      
      TYPE(meshType) :: m  ! mesh data structure
      TYPE(SuperElType) part(inp%maxNp),p
      TYPE(fluidType) :: fluid
      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position   

      INTEGER i,n,ierr
      REAL(8) a,b,c,e1,e2,rho,u(3)
      CHARACTER*64 patchName

      INTEGER ibound,sf,ef
      INTEGER iface, itri

      REAL(8) xp(3),xb(3),xtri(3,3),xcorr(3)

      
      REAL(8) r(3),v(3),o(3),x,phi,theta,psi,pi


      DO i = 1,m%nbound 

        !PRINT *, "bound(",i,")=",m%bound(i)%name 
        IF ( TRIM(m%bound(i)%name).EQ.TRIM(patchName) ) ibound = i

      END DO

      PRINT *, "Found patch at i =",i

      !Start and end face
      sf = m%bound(ibound)%startFace
      ef = sf + m%bound(ibound)%nFaces - 1

      !Loop for number of particles 
      DO i = 1,n

        p = part(i)

10      CONTINUE 

        !Random face
        CALL RANDOM_NUMBER(x)
        iface = sf + INT( x*(ef-sf) )
        !PRINT *, "iface=",iface

        !Set initial element from face owner
        p%element = m%f(iface)%owner

        !Random tri
        CALL RANDOM_NUMBER(x)
        itri = 1 + INT( x*(size(m%f(iface)%tf)-1) )

        !Random barycentrics
        CALL RANDOM_NUMBER(xb)
        xb = xb / sum(xb)

        xtri(1,:) = m%x( m%f(iface)%tf(itri)%con(1),: ) 
        xtri(2,:) = m%x( m%f(iface)%tf(itri)%con(2),: ) 
        xtri(3,:) = m%x( m%f(iface)%tf(itri)%con(3),: ) 

        xp = xb(1)*xtri(1,:) + xb(2)*xtri(2,:) + xb(3)*xtri(3,:)

        !Move just a little towards element center, to ensure it's found
        xcorr = 1.0E-03 * ( m%e( p%element )%xc - m%f(iface)%tf(itri)%center)
        p%r = xp + xcorr

        !Try if found, else go back, and try again with another one
        CALL GetFFFap(m,p,fluid,fap,-1,ierr)
        IF (ierr.NE.0) THEN 
          PRINT *, "ierr ne 0, trying again"
          GOTO 10
        ELSE
          !WRITE(*,'(A4,3F10.5)') "fap=",fap%vx
        END IF

        !WRITE(*,'(A4,3F10.5)') "x1=",xtri(1,:)
        !WRITE(*,'(A4,3F10.5)') "x2=",xtri(2,:)
        !WRITE(*,'(A4,3F10.5)') "x3=",xtri(3,:)
        !WRITE(*,'(A4,4F10.5)') "xb=",xb,sum(xb)
        !PRINT *, "p%el",p%element
        !WRITE(*,'(A4,3F10.5)') "p%r=",p%r

        !Set particle orientation
        v=0.0D0
        pi=4.0D0*ATAN(1.0D0)
        CALL RANDOM_NUMBER(x)
        theta=ACOS(2.0D0*x-1.0D0)  ! -pi .. pi
        CALL RANDOM_NUMBER(x)
        phi=2.0D0*pi*x ! ni fajn ce je nic
        CALL RANDOM_NUMBER(x)
        psi=2.0D0*pi*x
        o=0.0D0

        r = p%r
        CALL seInit2(p,a,b,c,e1,e2,rho,phi,theta,psi,r,v,o)
        !PRINT *, "after seinit"
        !PRINT *, "p%e=",p%element
     
        !If parsed nonzero, set parsed velocity
        IF ( u(1).NE.0.0D0.OR.u(2).NE.0.0D0.OR.u(3).NE.0.0D0 ) THEN
          p%v(1)=u(1)
          p%v(2)=u(2)
          p%v(3)=u(3)
        
        !If parsed zero, set to flow velocity
        ELSE
          CALL GetFFFap(m,p,fluid,fap,-1,ierr)
          p%v(1)=fap%vx
          p%v(2)=fap%vy
          p%v(3)=fap%vz
    
          p%vOld = p%v
    
          !WRITE(*,'(A4,3F10.5)') "p%r=",p%r
          !WRITE(*,'(A4,3F10.5)') "p%v=",p%v
        END IF
  
        p%pid=par_NextPID(cnt%LastPID)
        p%id=p%pid*100000000+cnt%cnp
        cnt%LastPID = p%pid
        
        IF (par_IsMine(p%pid)) THEN 
          cnt%cnp=cnt%cnp+1
          part(cnt%cnp)=p
        END IF

      END DO

      !CALL stopProgram(1)

      ierr = 0

      END 


! -----------------------------------------------------------------------------------------
      SUBROUTINE AddParticleKnownPosition(m,fluid,part,n,a,b,c,e1,e2,rho,box,ierr)
!
!     Add random particles
!      
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE mesh
      USE logFile
      USE counters
      USE ParalelEnvironment
      USE mFluid

      IMPLICIT NONE
      
      TYPE(meshType) :: m  ! mesh data structure
      TYPE(SuperElType) part(inp%maxNp),p
      TYPE(fluidType) :: fluid
      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position   


      INTEGER i,n,ierr,ie,ie2,itet,iface
      LOGICAL found

      REAL(8) dx,dxMin
      REAL(8) a,b,c,e1,e2,rho,box(3,2)
      REAL(8) r(3),v(3),o(3),x,phi,theta,psi,pi

      found = .FALSE.
      itet = 0
      ierr = 1

      p%r(1) = box(1,1)
      p%r(2) = box(1,2)
      p%r(3) = box(2,1)

      p%rOld = p%r

      PRINT *, "Starting distance calc"

      !Find closest element by distance to element center
      dxMin=9.9E+10;
      DO i = 1,m%nelem

        CALL dist2P(p%r,m%e(i)%xc,dx)

        IF (dx.LT.dxMin) THEN
            ie = i
            dxMin = dx
        END IF

      END DO

      PRINT *, "Closest element is: ",ie
      PRINT *, "Element center is: ", m%e(ie)%xc
      PRINT *, "Element nodes:"
      DO i = 1,size( m%e(ie)%nodeCon )
        WRITE(*,'(A4,3F14.8)') "x=",m%x( m%e(ie)%nodeCon(i),: )
      END DO


      !Check if in closest element
      found = findWithinElement(p%r,ie,itet)
      IF(found) GOTO 10

      ie2 = ie

      !Check if in neighbours
      DO i = 1,size( m%e(ie2)%faceCon )

        iface = m%e(ie2)%faceCon(i)

        IF (iface.LT.0) THEN
          ie = m%f(abs(iface))%owner
        ELSE 
          ie = m%f(abs(iface))%neighbour
        END IF

        IF (ie.GT.0) THEN

          found = findWithinElement(p%r,ie,itet)
          IF(found) GOTO 10

        END IF

      END DO

      PRINT *, "AddParticleKnownPosition :: Particle not found, return ierr = 1"
      p%active=.FALSE.
      GOTO 13

10    CONTINUE

      !Set particle element
      p%element = ie

      PRINT *, "AddParticleKnownPosition :: Found in element", p%element


!
!     Set particle orientation
!
      v=0.0D0
      pi=4.0D0*ATAN(1.0D0)

      CALL RANDOM_NUMBER(x)
      theta=ACOS(2.0D0*x-1.0D0)  ! -pi .. pi
      CALL RANDOM_NUMBER(x)
      phi=2.0D0*pi*x ! ni fajn ce je nic
      CALL RANDOM_NUMBER(x)
      psi=2.0D0*pi*x

      o=0.0D0

      r = p%r
      CALL seInit2(p,a,b,c,e1,e2,rho,phi,theta,psi,r,v,o)
      !PRINT *, "after seinit"
      !PRINT *, "p%e=",p%element
!
!     Set particle velocity to flow velocity
!     

      CALL GetFFFap(m,p,fluid,fap,-1,ierr)

      !PRINT *, "fap=",fap%vx,fap%vy,fap%vz
      p%v(1)=fap%vx
      p%v(2)=fap%vy
      p%v(3)=fap%vz

      p%vOld = p%v

      !WRITE(*,'(A4,3F10.5)') "p%r=",p%r
      !WRITE(*,'(A4,3F10.5)') "p%v=",p%v

      p%pid=par_NextPID(cnt%LastPID)
      p%id=p%pid*100000000+cnt%cnp
      cnt%LastPID = p%pid
      
      IF (par_IsMine(p%pid)) THEN 
        cnt%cnp=cnt%cnp+1
        part(cnt%cnp)=p
      END IF

      ierr = 0

13    CONTINUE
      ierr = 1


      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------


      LOGICAL FUNCTION findWithinElement(r,ie,itet)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: ie
        REAL(8), INTENT(IN) :: r(3)
        INTEGER, INTENT(OUT) :: itet

        INTEGER i,ntet,ierr
        REAL(8) fCenters(4,3), fNormals(4,3),res

        PRINT *, "checking element: ", ie

        ntet = size (m%e( ie )%tet)
        DO i = 1,ntet
  
          fCenters = m%e( ie )%tet(i)%centers
          fNormals = m%e( ie )%tet(i)%normals
  
          CALL PartInElement2(r,fCenters,fNormals,4,res,ierr)
  
          IF (ierr.EQ.0) THEN
            itet = i
            findWithinElement = .TRUE.
            RETURN
          END IF

        END DO

        findWithinElement = .FALSE.
    
      END FUNCTION findWithinElement


      END 


! -----------------------------------------------------------------------------------------
      SUBROUTINE ExportOnePartToLog(p)
!
!     Write particle info to log file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE logFile
      IMPLICIT NONE
      CHARACTER OneLine*255
      TYPE(SuperElType) p

      CALL logWrite("One particle data:")
      CALL logVecWrite("  - Position                  : ",p%r,3,'(A6,A,3(D11.5,1X))')
      CALL logVecWrite("  - Velocity                  : ",p%v,3,'(A6,A,3(D11.5,1X))')
      CALL logVecWrite("  - Angular velocity in PFR   : ",p%o,3,'(A6,A,3(D11.5,1X))')
      CALL logVecWrite("  - Euler angles              : ",p%ea,3,'(A6,A,3(D11.5,1X))')
      CALL logVecWrite("  - Euler parameters          : ",p%ep,4,'(A6,A,4(D11.5,1X))')
      CALL logRealWrite("  - a semi axis [m]           : ",p%a,'(A6,A,D11.5)')
      CALL logRealWrite("  - b semi axis [m]           : ",p%a,'(A6,A,D11.5)')
      CALL logRealWrite("  - c semi axis [m]           : ",p%a,'(A6,A,D11.5)')
      CALL logRealWrite("  - lambda = b/a              : ",p%lambda,'(A6,A,D11.5)')
      CALL logRealWrite("  - e1 super ellipsoid par.   : ",p%e1,'(A6,A,D11.5)')
      CALL logRealWrite("  - e2 super ellipsoid par.   : ",p%e2,'(A6,A,D11.5)')
      CALL logRealWrite("  - density  [kg/m^3]         : ",p%density,'(A6,A,D11.5)')
      CALL logRealWrite("  - volume [m^3]              : ",p%volume,'(A6,A,D11.5)')
      CALL logRealWrite("  - mass [kg]                 : ",p%mass,'(A6,A,D11.5)')
      CALL logRealWrite("  - A (density ratio)         : ",p%AA,'(A6,A,D11.5)')
      CALL logRealWrite("  - R (density ratio)         : ",p%RR,'(A6,A,D11.5)')
      CALL logVecWrite("  - Settling velocity [-]      : ",p%vs,3,'(A6,A,3(D11.5,1X))')
      CALL logRealWrite("  - tau_particle [s]          : ",p%tau,'(A6,A,D11.5)')
      CALL logRealWrite("  - Stokes number             : ",p%St,'(A6,A,D11.5)')
      CALL logRealWrite("  - Elliptic tau_particle     : ",p%taue,'(A6,A,D11.5)')
      CALL logRealWrite("  - Elliptic St. number       : ",p%Ste,'(A6,A,D11.5)')
      CALL logVecWrite("  - Resistance tensor in PFR  : ",p%ResTprime,3,'(A6,A,3(D11.5,1X))')
      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE seBaseInit(se)
!
!     Initiallize super ellipsoid partice
!
! -----------------------------------------------------------------------------------------
      USE inpFile
      USE logFile
      USE SuperE      
      IMPLICIT NONE

      TYPE(SuperElType) se
      REAL(8) f
      INTEGER i


      se%lambda=se%c/se%a

      se%AA = se%density / (se%density + 0.5D0 * inp%f_rho)
      se%RR = inp%f_rho  / (0.5D0 * inp%f_rho + se%density)

      CALL seSetVolume(se)
      CALL seMOI(se)
      se%mass=se%density*se%volume

!
!     Particle response time
!
      se%tau = 2.0D0 * se%density * se%a * se%a / ( 9.0D0 * inp%f_rho * inp%f_mu )
      se%St = se%tau / (inp%f_L/inp%f_u0)
!
!     Elliptic particle response time
!      
      IF (se%lambda.LT.1.0D0) THEN
        CALL logWrite("ERROR :: seInit :: Particle aspect ratio < 1 !")   
        CALL StopProgram(1)   
      ELSE IF (se%lambda.EQ.1.0D0) THEN
        f = 1.0D0
      ELSE
        f = se%lambda*Log(se%lambda+sqrt(se%lambda**2-1.0D0)) / sqrt(se%lambda**2-1.0D0)
      END IF
      se%Ste=se%st*f
      se%taue=se%tau*f
!
!     Settling velocity
!
      DO i=1,3
        se%vs(i) = se%tau * inp%g(i) * ( 1.0D0 - inp%f_rho / se%density ) / inp%f_u0
      END DO
!
!     Calculate rotation matrix
!
      CALL seCalRotationMatrix(se)
!
!     Particle orientation in GRF
!
      CALL seOrientationAxis(se)
!
!     Calculate size of shpere enclosing the super ellipsoid
!
      CALL seCalVOS(se)
!
!     Calculate rotation dynamics eq. factors (rdfx)
!
      CALL CalRotDynFactor(se%rdfx,se%rdfy,se%rdfz,se%a,se%lambda, &
                                 inp%f_mu,inp%f_rho,inp%f_L,se%density,inp%f_u0)
!
!     Particle resistance tensor (prime)
!
      CALL CalResTensorPrime(se%lambda,se%ResTprime)
!
!     Get euler angels from parameters
!
      CALL seEulerParameters2AnglesRotMat(se)

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE seInit(se,a,b,c,e1,e2,rho,phi,theta,psi,r,v,o)
!
!     Initiallize super ellipsoid partice
!
! -----------------------------------------------------------------------------------------
      USE inpFile
      USE logFile
      USE SuperE      
      IMPLICIT NONE

      TYPE(SuperElType) se
      REAL(8) a,b,c,e1,e2,rho,phi,theta,psi,f
      REAL(8) r(3),v(3),o(3)
      INTEGER i

      se%a=a
      se%b=b
      se%c=c
      se%e1=e1
      se%e2=e2
      se%density=rho

      se%ea(1)=psi
      se%ea(2)=theta
      se%ea(3)=phi

      CALL seEulerAngles2Parameters(se)
      se%r=r
      se%v=v
      se%o=o
      se%element=1  ! dont know where it is yet
      se%active=.TRUE.
!
!     Do basic init from data in se
!
      CALL seBaseInit(se)

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE seInit2(se,a,b,c,e1,e2,rho,phi,theta,psi,r,v,o)
!
!     Initiallize super ellipsoid partice
!
! -----------------------------------------------------------------------------------------
      USE inpFile
      USE logFile
      USE SuperE      
      IMPLICIT NONE

      TYPE(SuperElType) se
      REAL(8) a,b,c,e1,e2,rho,phi,theta,psi,f
      REAL(8) r(3),v(3),o(3)
      INTEGER i

      se%a=a
      se%b=b
      se%c=c
      se%e1=e1
      se%e2=e2
      se%density=rho

      se%ea(1)=psi
      se%ea(2)=theta
      se%ea(3)=phi

      CALL seEulerAngles2Parameters(se)
      se%r=r
      se%v=v
      se%o=o
      se%active=.TRUE.
!
!     Do basic init from data in se
!
      CALL seBaseInit(se)

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE CalResTensorPrime(lambda,ResTprime)
!
!     Calculate resistance tensor in particle frame of reference
!
! -----------------------------------------------------------------------------------------
      USE logFile
      IMPLICIT NONE
      REAL(8) lambda,ResTprime(3),l2

      ResTprime=0.0D0
      IF (lambda.LT.1.0D0) THEN
        CALL logWrite("ERROR :: CalResTensorPrime :: Particle aspect ratio < 1 !")   
        CALL StopProgram(1)
      ELSE IF (lambda.EQ.1.0D0) THEN
	    ResTprime(1)=6.0D0
	    ResTprime(2)=6.0D0
	    ResTprime(3)=6.0D0
      ELSE
        l2=Lambda**2 - 1.0D0
    	  ResTprime(1)=(16.0D0*l2**(1.5D0))/( (2*l2 - 1.0D0)*Log(Lambda + Sqrt(l2)) + Lambda*Sqrt(l2) )
        ResTprime(2)=ResTprime(1)
        ResTprime(3)=( 8.0D0*l2**(1.5D0))/( (2*l2 + 1.0D0)*Log(Lambda + Sqrt(l2)) - Lambda*Sqrt(l2) )
      END IF

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE ParticleListExportCSV(part)
!
!     Write to csv file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE logFile
      USE mesh
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) part(inp%maxNp)
      TYPE(localMeshType)  :: sphere
      REAL(8) time,r(3)
      INTEGER np,i,j,k,lun
      CHARACTER(1000) vrstica
 

      lun = 12

!
!     Create subdirectory for Paraview files
!
!      CALL create_directory( inp%ParaviewResultsFolder )
!
!     Determine file name
!
      WRITE(vrstica,'(A,A,A,A,I0)') TRIM(inp%ParaviewResultsFolder),"/",TRIM(lg%IDname),".para-list.csv.",cnt%iTime
      OPEN (lun,FILE=TRIM(vrstica),ERR=10,STATUS='UNKNOWN')
      WRITE (lun,'(A)') 'x,y,z,vx,vy,vz,ox,oy,oz,axisX,axisY,axisZ,ea1,ea2,ea3,ep1,ep2,ep3,ep4,ierr,element'
!
!     Loop over particles to export
!
      DO i=1,cnt%cnp
!
!       Write to file
!
        WRITE(vrstica,'(20(G20.10,A))') part(i)%r(1),",",part(i)%r(2),",",part(i)%r(3),"," &
     &                                 ,part(i)%v(1),",",part(i)%v(2),",",part(i)%v(3),"," &
     &                                 ,part(i)%o(1),",",part(i)%o(2),",",part(i)%o(3),","      &
     &                                 ,part(i)%axis(1),",",part(i)%axis(2),",",part(i)%axis(3),","      &
     &                                 ,part(i)%ea(1),",",part(i)%ea(2),",",part(i)%ea(3),","      &
                                      ,part(i)%ep(1),",",part(i)%ep(2),",",part(i)%ep(3),",",part(i)%ep(4),",",part(i)%ierr,","   
        IF (part(i)%active) THEN
             WRITE(vrstica,'(A,I0)') TRIM(vrstica),part(i)%element
        ELSE
             WRITE(vrstica,'(A,I0)') TRIM(vrstica),-part(i)%element
            
        END IF
        CALL sqblnk(lun,vrstica)
      END DO

      CLOSE (lun)

      RETURN

10    CONTINUE ! error when opening input file
      CALL logWrite ("ERROR :: ParticleListExportCSV :: Could not open file for export!")
      CALL StopProgram(1)

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE ParticleListExportVTK(part)
!
!     Write to VTK file
!
! -----------------------------------------------------------------------------------------
      USE superE
      USE inpFile
      USE logFile
      USE mesh
      USE counters
      IMPLICIT NONE

      TYPE(SuperElType) part(inp%maxNp)
      TYPE(localMeshType)  :: sphere
      REAL(8) time,r(3)
      INTEGER np,i,j,k,lun
      CHARACTER(1000) vrstica
 

      lun = 12
!
!     Create subdirectory for Paraview files
!
!      CALL create_directory( inp%ParaviewResultsFolder )
!
!     Determine file name
!
      WRITE(vrstica,'(A,A,A,A,I0)') TRIM(inp%ParaviewResultsFolder),"/",TRIM(lg%IDname),".para-list.vtk.",cnt%iTime
      OPEN (lun,FILE=TRIM(vrstica),ERR=10,STATUS='UNKNOWN')
      WRITE (lun,'(A)') '# vtk DataFile Version 2.0'
      WRITE (lun,'(A)') 'superV results file'
      WRITE (lun,'(A)') 'ASCII'

      WRITE (lun,'(A)') 'DATASET UNSTRUCTURED_GRID'

!
!     Locations
!
      WRITE (lun,'(A,I0,A)') 'POINTS ',cnt%cnp,' float'
      DO i=1,cnt%cnp
            WRITE(vrstica,'(3(G20.10))') part(i)%r(1),part(i)%r(2),part(i)%r(3)
            CALL sqblnk(lun,vrstica)
      END DO
!
!     Cell types = 1 = vertex
!
      WRITE (lun,'(A,I0,A)') 'CELL_TYPES ',cnt%cnp
      DO i=1,cnt%cnp
            WRITE(vrstica,'(A)') "1"
            CALL sqblnk(lun,vrstica)
      END DO

      WRITE (lun,'(A,I0)') 'POINT_DATA ',cnt%cnp

      WRITE (lun,'(A)') 'VECTORS Velocity float'
      DO i=1,cnt%cnp
            WRITE(vrstica,'(3(G20.10))') part(i)%v(1),part(i)%v(2),part(i)%v(3)
            CALL sqblnk(lun,vrstica)
      END DO

      WRITE (lun,'(A)') 'VECTORS AngularVelPFR float'
      DO i=1,cnt%cnp
            WRITE(vrstica,'(3(G20.10))') part(i)%o(1),part(i)%o(2),part(i)%o(3)
            CALL sqblnk(lun,vrstica)
      END DO

      WRITE (lun,'(A)') 'VECTORS Axis float'
      DO i=1,cnt%cnp
            WRITE(vrstica,'(3(G20.10))') part(i)%axis(1),part(i)%axis(2),part(i)%axis(3)
            CALL sqblnk(lun,vrstica)
      END DO

      WRITE (lun,'(A)') 'VECTORS EulerAngles float'
      DO i=1,cnt%cnp
            WRITE(vrstica,'(3(G20.10))') part(i)%ea(1),part(i)%ea(2),part(i)%ea(3)
            CALL sqblnk(lun,vrstica)
      END DO

      WRITE (lun,'(A)') 'SCALARS element integer 1'
      WRITE (lun,'(A)') 'LOOKUP_TABLE default'
      DO i=1,cnt%cnp
        IF (part(i)%active) THEN
             WRITE(vrstica,'(I0)') part(i)%element
        ELSE
             WRITE(vrstica,'(I0)') -part(i)%element
            
        END IF
        CALL sqblnk(lun,vrstica)
      END DO

      WRITE (lun,'(A)') 'SCALARS active integer 1'
      WRITE (lun,'(A)') 'LOOKUP_TABLE default'
      DO i=1,cnt%cnp
        IF (part(i)%active) THEN
             WRITE(vrstica,'(I0)') 1
        ELSE
             WRITE(vrstica,'(I0)') 0
            
        END IF
        CALL sqblnk(lun,vrstica)
      END DO

      CLOSE (lun)

      RETURN

10    CONTINUE ! error when opening input file
      CALL logWrite ("ERROR :: ParticleListExportVKT :: Could not open file for export!")
      CALL StopProgram(1)

      END

