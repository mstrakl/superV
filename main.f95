!
!     SuperV - superellipsoid particle tracking
!
      PROGRAM SSS
!      
!     Written by
!
!           - Jure Ravnik, jure.ravnik@um.si
!
!           University of Maribor
!           Faculty of Mechanical Engineering
!           Smetanova 17, SI-2000, Maribor, Slovenia
!
!     Web page
!
!          http://jure.ravnik.si
!
!
!     Description
!
!     Vreteno is a software simulation tool capable of simualting movement of
!     particles in fluid flows under the influence of hydrodynamic and magnetic forces
!
!
!     Licence
!
!     This software may be freely used provided that any results, that were obtained by the
!     use of this software or parts of it, cite the following papers:
!
!     RAVNIK, Jure, HRIBERŠEK, Matjaž. High gradient magnetic particle separation in viscous flows by 3D BEM.
!     Comput. mech., Online First, 24 May 2012, doi: 10.1007/s00466-012-0729-3.
!     http://dx.doi.org/10.1007/s00466-012-0729-3
!
!     RAVNIK, Jure, ŠKERGET, Leopold, HRIBERŠEK, Matjaž, ŽUNIČ, Zoran.
!     Numerical simulation of dilute particle laden flows by wavelet BEM-FEM.
!     Comput. methods appl. mech. eng.. [Print ed.], Jan. 2008, vol. 197, iss. 6/8, pp. 789-805.
!     http://dx.doi.org/10.1016/j.cma.2007.09.007.
!
!     where the theory behind this software is described.
!
!
!     Requirements
!
!
!     Makefile is written for Inter Fortran Compiler, but the code should work
!     with any FORTRAN 90 compiler. When linking, MPI library must be likned,
!     version 1.2.7 is used by the author. See Makefile in source folder for
!     further information.
!
!
!     Disclaimer of Warranty.
!
!     THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY
!     APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT
!     HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM "AS IS" WITHOUT WARRANTY
!     OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,
!     THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
!     PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM
!     IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF
!     ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
!
!     Limitation of Liability.
!
!     IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
!     WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MODIFIES AND/OR CONVEYS
!     THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY
!     GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE
!     USE OR INABILITY TO USE THE PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF
!     DATA OR DATA BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD
!     PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS),
!     EVEN IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
!     SUCH DAMAGES.
!
      USE SuperE
      USE logFile
      USE inpFile
      USE counters
      USE mesh
      USE mFluid
      USE cpuTime
      USE ParalelEnvironment
      IMPLICIT NONE

      TYPE(SuperElType), ALLOCATABLE :: part(:)
      TYPE(meshType) :: FoamMesh
      TYPE(fluidType) :: fluid
      TYPE(fluidAtPartType) :: fap ! fluid flow fields at particle position

      INTEGER i,j,iT,ierr
 
      CALL cpInit()
      CALL cpStart(cput%meas(1))
! 
!     Version data 
!
      lg%IDname='superV'
      lg%IDversion='0.7'
      lg%IDdate='Nov 2020'
      CALL RANDOM_SEED      
!
!     Init paralel environment
!
      CALL par_init()    
!
!     Open log file + init
!
      CALL logInit()
!
!     Read Input File
!
      CALL inpReadInputFile()     
!
!     Read mesh file
!
      CALL GetVTKmesh(FoamMesh,inp%VTKflowResults,inp%FoamMeshDir)
!
!     Allocate memory for fluid flow field
!
      PRINT *, "Allocating fluid fields"
      CALL AllocateFluidFileds(FoamMesh,fluid)
!
!     Init time
!
      CALL cntInit()
!
!     Get flow field
!
      PRINT *, "Reading flow field"
      CALL GetFlowField(FoamMesh,fluid) ! PAZI : fluid%dvxdt ni izračunan = 0 !
!
!     Allocate memory for particles
!
      IF (inp%maxNp.EQ.0) THEN
        PRINT *, "Start EstimateNumberOfParticles"
        CALL EstimateNumberOfParticles()
        CALL logIntWrite("Total number of particles : ",inp%maxNp)
      ELSE
        CALL logIntWrite("Memory reserved for particles : ",inp%maxNp)
      END IF
      ALLOCATE (part(inp%maxNp))
!
!     Add particles to model
!
      CALL logWrite("Adding initial particles")
      PRINT *, "Start AddParticles"
      CALL AddParticles(FoamMesh,part,fluid)

      PRINT *, "Initial number of particles is:", cnt%cnp

      CALL BreakPoint()

      IF (cnt%cnp.LT.1) THEN
        PRINT *, "Stopping after AddParticles, as there are no particles"
        CALL stopProgram(1)
      END IF
!
!     Export 1 particle to log file
!      
      CALL ExportOnePartToLog(part(1))
!
!     Export initial particle distribution to files
!
      CALL ExportToBIN(part)
!
!     Time loop
!
      CALL logWrite ("Start time loop")
      CALL cpStart(cput%meas(2))            
      DO iT=1,inp%nTimeSteps
        CALL logPercent(iT,inp%nTimeSteps,10)
        cnt%rTime = cnt%rTime + inp%TimeStep
        cnt%iTime = cnt%iTime + 1

        PRINT *,"------------"
        PRINT *,"iT=",iT," Lost=",cnt%lost


!
!       For all active particles
!           
        DO i=1,cnt%cnp


            !PRINT *, "---------------------------------------"
            !WRITE(*,'(A10,I6)') "p%elem=",part(i)%element
            !WRITE(*,'(A10,3F12.6)') "fap%v=",fap%vx,fap%vy,fap%vz
            !WRITE(*,'(A10,3F12.6)') "p%r=",part(i)%r
            !WRITE(*,'(A10,3F12.6)') "p%rOld=",part(i)%rOld
            !WRITE(*,'(A10,3F12.6)') "p%v=",part(i)%v
            !WRITE(*,'(A10,3F12.6)') "p%vOld=",part(i)%vOld

          !PRINT *, "start interpolation testing"
          !CALL interpolationTesting(VTKmesh,part(i),fluid,fap,iT,ierr)
          !PRINT *, "end interpolation testing"
          !CALL stopProgram(1)

          !PRINT *, "start grad testing"
          !CALL gradientTesting(VTKmesh,part(i),fluid,fap,iT,ierr)
          !PRINT *, "end grad testing"


          IF (part(i)%active) THEN
!
!           Get flow field data at the position of the particle
!        
            CALL GetFFFap(FoamMesh,part(i),fluid,fap,iT,ierr)
!
!           Check if particle is outside the mesh
!
            IF (ierr.NE.0) THEN
!
!             Particle hit boundary
!             
              IF (ierr.EQ.2) THEN    
                part(i)%ierr = ierr
                part(i)%active=.FALSE.
              END IF

              IF (ierr.EQ.1) THEN
                  !IF(part(i)%nnf.GE.10) THEN
                  part(i)%ierr = ierr
                  part(i)%active=.FALSE.
                  cnt%lost = cnt%lost + 1

                  PRINT *, "------------------------------"
                  PRINT *, "Particle lost at iT=",iT
                  PRINT *, "Part ID=",part(i)%id
                  PRINT *, "Element ID=",part(i)%element
                  PRINT *, "------------------------------"

                  CALL WriteVTKTriangulatedElement('diag/trielement',iT,FoamMesh,part(i)%element)
                  CALL WriteVTKPoint('diag/point-cc',iT,FoamMesh%e(part(i)%element)%xc)
                  CALL WriteVTKPoint('diag/point-a',iT,part(i)%rOld)
                  CALL WriteVTKPoint('diag/point-b',iT,part(i)%r)

                  !END IF
              END IF

            ELSE
!
!             Move particle
!
              CALL  MoveParticle(part(i),fap,cput%meas(7))
              !PRINT *, "---------------------------------------"
              !WRITE(*,'(A10,I6)') "p%elem=",part(i)%element
              !WRITE(*,'(A10,3F12.6)') "fap%v=",fap%vx,fap%vy,fap%vz
              !WRITE(*,'(A10,3F12.6)') "p%r=",part(i)%r
              !WRITE(*,'(A10,3F12.6)') "p%rOld=",part(i)%rOld
              !WRITE(*,'(A10,3F12.6)') "p%v=",part(i)%v
              !WRITE(*,'(A10,3F12.6)') "p%vOld=",part(i)%vOld

              IF (inp%Periodic.GT.0) CALL PeriodicBC(part(i))
            END IF
          END IF
        END DO  

!
!       Add particles to model
!  
        CALL AddParticles(FoamMesh,part,fluid)
!
!       Export particles to BIN files
!
        IF (cntExport(inp%ResultExportFreq)) CALL ExportToBIN(part)
!  
!       Stop for debugging
!
        !IF (iT.EQ.5) THEN
        !  PRINT *, "stopping main loop, after timestep:",iT
        !  CALL stopProgram(1)
        !END IF

      END DO

12    continue
      CALL logWrite ("End time loop")
      CALL cpStop(cput%meas(2))      
      CALL cpStop(cput%meas(1))

!     
!     Blesavi printi
!     
      IF (par_AmIRankZero()) THEN

      PRINT *,"--------------------------------------"
      PRINT *, "Interpolate Particle Search:"
      PRINT *, "found correctly=",cnt%inthowfound(1)
      PRINT *, "not found=",cnt%inthowfound(2)
      PRINT *, "found multi=",cnt%inthowfound(3)

      PRINT *,"--------------------------------------"
      PRINT *, "Particles lost=",cnt%lost

!      PRINT *,"--------------------------------------"
!      PRINT *,"Particle element history:"
!
!      DO i = 1,cnt%cnp
!            j=1
!            DO, WHILE ( part(i)%cellHist(j,1).NE.-1 )
!            PRINT *, "p(",i,")%ch(",j,")=",part(i)%cellHist(j,:)
!            j = j + 1
!            END DO
!      END DO


      ENDIF
!
!     Export to ASCII
!
      IF (par_AmIRankZero()) CALL ExportResultsToASCII(part,FoamMesh,fluid)
!
!     Write CPU time measurements to log file
!
      CALL cpWriteToLog()
!
!     Close files and stop program
!
      CALL StopProgram(0)


      END
