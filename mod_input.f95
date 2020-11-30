!
!
!
      MODULE inpFile
!
! -----------------------------------------------------------------------------------------
!
      TYPE InpType
        INTEGER lun
        CHARACTER(100) FileName
        INTEGER  maxNp   ! maximal number of particles for memory allocation
        INTEGER ResultExportFreq ! when to write binary files
        INTEGER TecplotMeshExport ! switch that togles Tecplot meshed particle export
        INTEGER ParaviewMeshExport ! switch that togles Tecplot meshed particle export
        INTEGER ParaviewListExport ! switch paraview list CSV export
        INTEGER ParaviewListVTKExport ! switch paraview list VTK export
        INTEGER WriteVTKflowResults ! switch to write flow results
        INTEGER AnalFlowField ! is the flow field analytically computed
        INTEGER Periodic,PeriodicAxis
        REAL(8) PeriodicFrom,PeriodicTo
        INTEGER MoveModel ! =0 massless, =1 with mass
        CHARACTER(50) ParaviewResultsFolder ! where to store paraview files (there is a lot of them)
        CHARACTER(50) BINResultsFolder ! where to store bin files (there is a lot of them)
        INTEGER iTME ! first time counter (do not edit)
        CHARACTER(100) VTKflowResults ! path and filename of flow results or path to OpenFOAM folder
        CHARACTER(100) FoamMeshDir ! path to foam mesh dir, to read supporting files
        INTEGER nTimeSteps ! total number of time steps to perform
        REAL(8) TimeStep  ! time step size
        REAL(8) f_L        ! characteristic dimension for flow ( f_Re = f_l * f_u0 / f_mu )
        REAL(8) f_u0       ! characteristic velocity for flow ( f_Re = f_l * f_u0 / f_mu )
        REAL(8) f_mu       ! kinematic viscosity  ( f_Re = f_l * f_u0 / f_mu )
        REAL(8) f_rho      ! fluid density 
        REAL(8) f_Re       ! flow Reynolds number ( f_Re = f_l * f_u0 / f_mu )
        INTEGER nPartStat  ! number of particle statements in input file
        CHARACTER(255), ALLOCATABLE :: PartStat(:) ! part statements (nPartStat)
        REAL(8) g0,g(3)    ! gravity vector 
        INTEGER fm_Gravity,fm_StokesDrag,fm_EllipticDrag,fm_AmPc

        !Mitja added
        INTEGER PurgeOldRes(4)
        INTEGER maxPartNNF
        INTEGER polyInterpol !poly interpolation method

      END TYPE

      TYPE(InpType) :: inp

      CONTAINS

! -----------------------------------------------------------------------------------------
      SUBROUTINE inpReadInputFile()
!
!     Read Input File
!
! -----------------------------------------------------------------------------------------
      USE logFile
      USE String_Utility
      USE ParalelEnvironment
      IMPLICIT NONE
      CHARACTER KeyWord*64,OneLine*255,dummy*64
      INTEGER i

      inp%lun = 12
      WRITE (inp%FileName,'(A,A)') TRIM(lg%IDname),".inp"

      OPEN (inp%lun,FILE=inp%FileName,ERR=10,STATUS='OLD')

      CALL logWrite ("MESSAGE :: ReadInputFile :: Reading : "//TRIM(inp%FileName))

!
!***    Set up default settings:
!
      inp%maxNp=0
      inp%ResultExportFreq=0
      inp%TecplotMeshExport=0
      inp%ParaviewMeshExport=0
      inp%ParaviewListExport=0
      inp%WriteVTKflowResults=0
      inp%ParaviewListVTKExport=0
      inp%MoveModel=0
      inp%iTME = 1
      inp%TimeStep=1.0D0
      inp%ParaviewResultsFolder="parares"
      inp%BINResultsFolder="binary"
      inp%VTKflowResults="NODATA"
      inp%FoamMeshDir="NODATA"
      inp%AnalFlowField=0
      inp%f_L = 1.0D0       ! characteristic dimension for flow ( f_Re = f_l * f_u0 / f_mu )
      inp%f_u0= 1.0D0       ! characteristic velocity for flow ( f_Re = f_l * f_u0 / f_mu )
      inp%f_mu= 1.0D0       ! kinematic viscosity  ( f_Re = f_l * f_u0 / f_mu )
      inp%f_rho= 1.0D0      ! fluid density 
      inp%f_Re= 1.0D0       ! flow Reynolds number ( f_Re = f_l * f_u0 / f_mu )
      inp%nPartStat = 0
      inp%fm_Gravity=0
      inp%fm_StokesDrag=0
      inp%fm_EllipticDrag=0
      inp%fm_AmPc=0
      inp%Periodic=0
      inp%PeriodicAxis=1
      inp%PeriodicFrom=0.0D0
      inp%PeriodicTo=0.0D0

      !Mitja added
      inp%PurgeOldRes=0
      inp%maxPartNNF=10
      inp%polyInterpol=1

!
!     Read
!
      CALL rOneTL(inp%lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        READ (OneLine,*) KeyWord
!
!       Reading keywords
!
        IF (StrLowCase(TRIM(KeyWord)).EQ.'maxnumberofparticles') THEN
          READ(OneLine,*) dummy,inp%maxNp
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'tecplotmeshexport') THEN
          READ(OneLine,*) dummy,inp%TecplotMeshExport
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'resultexportfreq') THEN
          READ(OneLine,*) dummy,inp%ResultExportFreq          
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'paraviewresultsfolder') THEN
          READ(OneLine,*) dummy,inp%ParaviewResultsFolder
          ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'binresultsfolder') THEN
          READ(OneLine,*) dummy,inp%BINResultsFolder          
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'paraviewmeshexport') THEN
          READ(OneLine,*) dummy,inp%ParaviewMeshExport
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'paraviewlistexport') THEN
          READ(OneLine,*) dummy,inp%ParaviewListExport
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'paraviewlistvtkexport') THEN
          READ(OneLine,*) dummy,inp%ParaviewListVTKExport
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'writevtkflowresults') THEN
          READ(OneLine,*) dummy,inp%WriteVTKflowResults
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'analflowfield') THEN
          READ(OneLine,*) dummy,inp%AnalFlowField     
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'vtkflowresults') THEN
          READ(OneLine,'(A15,A)') dummy,inp%VTKflowResults
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'foammeshdir') THEN
          READ(OneLine,'(A12,A)') dummy,inp%FoamMeshDir
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'timestep') THEN
          READ(OneLine,*) dummy,inp%TimeStep
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'periodic') THEN
          READ(OneLine,*) dummy,inp%Periodic,inp%PeriodicAxis,inp%PeriodicFrom,inp%PeriodicTo
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'numberoftimesteps') THEN
          READ(OneLine,*) dummy,inp%nTimeSteps
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'movemodel') THEN
          READ(OneLine,*) dummy,inp%MoveModel
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'forces') THEN
          READ(OneLine,*) dummy,inp%fm_Gravity,inp%fm_StokesDrag,inp%fm_EllipticDrag,inp%fm_AmPc
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'fluidflow') THEN
          READ(OneLine,*) dummy,inp%f_L,inp%f_u0,inp%f_mu,inp%f_rho
          inp%f_Re=inp%f_L*inp%f_u0/inp%f_mu
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'gravity') THEN
          READ(OneLine,*) dummy,inp%g0,inp%g(1),inp%g(2),inp%g(3)
          inp%g(1)=inp%g0*inp%g(1)
          inp%g(2)=inp%g0*inp%g(2)
          inp%g(3)=inp%g0*inp%g(3)
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'particles') THEN
          inp%nPartStat = inp%nPartStat + 1
!
!         Mitja added 
!
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'purgeoldres') THEN
          READ(OneLine,*) dummy,inp%PurgeOldRes
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'maxpartnnf') THEN
          READ(OneLine,*) dummy,inp%maxPartNNF
        ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'polyinterpol') THEN
          READ(OneLine,*) dummy,inp%polyinterpol
        ELSE
          CALL logWrite ("WARNING :: ReadInputFile :: Unknown keyword : "//TRIM(KeyWord))
        END IF
        CALL rOneTL(inp%lun,OneLine)
      END DO

!
!     Purge old results files
!
      IF (par_AmIRankZero()) THEN
        IF (inp%PurgeOldRes(1).EQ.1) &
        CALL SYSTEM ("if [ -f *.sVbin ]; then rm -r *.sVbin; fi")

        IF (inp%PurgeOldRes(2).EQ.1) &
        CALL SYSTEM ("if [ -f *.vtk.nec ]; then rm -r *.vtk.nec; fi")

        IF (inp%PurgeOldRes(3).EQ.1) &
        CALL SYSTEM ("if [ -d 'parares' ]; then rm -r parares; fi")

        IF (inp%PurgeOldRes(4).EQ.1) &
        CALL SYSTEM ("if [ -d 'binary' ]; then rm -r binary; fi")

        !PRINT *, "purge=",inp%purgeOldRes
        !PRINT *, "deleting"
      END IF

!     Read particle statements
      ALLOCATE (inp%PartStat(inp%nPartStat))
      inp%nPartStat = 0
      REWIND(inp%lun)
      CALL rOneTL(inp%lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        READ (OneLine,*) KeyWord
         IF (StrLowCase(TRIM(KeyWord)).EQ.'particles') THEN
          inp%nPartStat = inp%nPartStat + 1
          READ(OneLine,'(A10,A)') dummy,inp%PartStat(inp%nPartStat)
        END IF
        CALL rOneTL(inp%lun,OneLine)
      END DO

      CLOSE (inp%lun)


      CALL logWrite("Characteristic fluid & flow data: ")
      WRITE (OneLine,'(A,E10.5,A)') "  - Characteristic dimension : ",inp%f_L," m."
      CALL logWrite(TRIM(OneLine))
      WRITE (OneLine,'(A,E10.5,A)') "  - Characteristic velocity  : ",inp%f_u0," m/s."
      CALL logWrite(TRIM(OneLine))
      WRITE (OneLine,'(A,E10.5,A)') "  - Fluid kinematic viscosity: ",inp%f_mu," m^2/s."
      CALL logWrite(TRIM(OneLine))
      WRITE (OneLine,'(A,E10.5,A)') "  - Fluid dynamic viscosity  : ",inp%f_mu*inp%f_rho," Pa*s."
      CALL logWrite(TRIM(OneLine))
      WRITE (OneLine,'(A,E10.5,A)') "  - Fluid Density            : ",inp%f_rho," kg/m^3."
      CALL logWrite(TRIM(OneLine))
      WRITE (OneLine,'(A,E10.5,A)') "  - Fluid Reynolds number    : ",inp%f_Re,"."
      CALL logWrite(TRIM(OneLine))
      CALL logVecWrite             ("  - Gravity                  : ",inp%g,3,'(A6,A,3(D11.5,1X))')

      CALL logWrite("Particle statements: ")
      DO i=1,inp%nPartStat
        WRITE (OneLine,'(A,A)') "  - ",TRIM(inp%PartStat(i))
        CALL logWrite(TRIM(OneLine))
      END DO
      
      CALL logWrite("Forces considered: ")
      IF (inp%fm_Gravity.GT.0)      CALL logWrite("  - Gravity")
      IF (inp%fm_StokesDrag.GT.0)   CALL logWrite("  - Stokes drag")
      IF (inp%fm_EllipticDrag.GT.0) CALL logWrite("  - Elliptic drag")
      IF (inp%fm_AmPc.GT.0)         CALL logWrite("  - Added mass & pressure correction")

!
!     Create subdirectory for Paraview & bin files
!
      CALL create_directory( inp%ParaviewResultsFolder )      
      CALL create_directory( inp%BINResultsFolder )      
      

      RETURN
      
10    continue ! error when opening input file
      CALL logWrite ("ERROR :: ReadInputFile :: Could not open : "//TRIM(inp%FileName))
      CALL StopProgram(1)
      
      END SUBROUTINE



      END MODULE
