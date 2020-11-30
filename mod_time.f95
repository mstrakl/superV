!
!
!
      MODULE cpuTime
!
! -----------------------------------------------------------------------------------------
!
      IMPLICIT NONE

      TYPE cpuTimeType
        CHARACTER(100), ALLOCATABLE :: desc(:)
        REAL(8), ALLOCATABLE :: meas(:)
        INTEGER n 
      END TYPE

      TYPE(cpuTimeType) :: cput

      CONTAINS

! -----------------------------------------------------------------------------------------
      SUBROUTINE cpInit()
!
!     Init time measurement
!
! -----------------------------------------------------------------------------------------

      INTEGER i
      
      cput%n=8
      ALLOCATE (cput%desc(cput%n))
      ALLOCATE (cput%meas(cput%n))

      cput%desc(1) = "TOTAL"
      cput%desc(2) = "time loop"
      cput%desc(3) = "GetFFFap"
      cput%desc(4) = "FindParticleInMesh -"
      cput%desc(5) = "GetXiEtaZeta -"
      cput%desc(6) = "InterpolateFFInElement -"
      cput%desc(7) = "MoveParticle"
      cput%desc(8) = "SetGrad"



      DO i=1,cput%n
        cput%meas(i)=0.0D0
      END DO


      END SUBROUTINE



! -----------------------------------------------------------------------------------------
      SUBROUTINE cpWriteToLog()
!
!     Write timing report to log file
!
! -----------------------------------------------------------------------------------------
      USE logFile
      INTEGER i

      CALL logWrite("")
      CALL logWrite("--------- CPU time measurements ---------")
      DO i=1,cput%n
        CALL logRealWrite(cput%desc(i),cput%meas(i),'(A6,A26,F15.3)')
      END DO
      CALL logWrite("")

      END SUBROUTINE



! -----------------------------------------------------------------------------------------
      SUBROUTINE cpStart(c)
!
!     Set time to zero
!
! -----------------------------------------------------------------------------------------
      REAL(8) c
      c = cpTime(0.0D0)

      END SUBROUTINE


! -----------------------------------------------------------------------------------------
      SUBROUTINE cpStop(c)
!
!     Set time to time elapsed
!
! -----------------------------------------------------------------------------------------
      REAL(8) c
      c = cpTime(c)

      END SUBROUTINE

!----------------------------------------------------------------------C

      REAL(8) FUNCTION cpTime (rvar)

!----------------------------------------------------------------------C
! **********************************************************************
! **                                                                  **
! ** Returns the CPU TIME in seconds                                  **
! **             --  ----                                             **
! **********************************************************************
!
      REAL(8) rvar,tmp
      
      CALL CPU_TIME(tmp)
      cptime = tmp - rvar

      END FUNCTION


      END MODULE      