!
!
!
      MODULE counters
!
! -----------------------------------------------------------------------------------------
!
      TYPE CntType
        INTEGER iTime    ! time step number
        REAL(8) rTime    ! real time
        INTEGER LastPID  ! ID of processor that got the last particle 
        INTEGER cnp      ! current number of particles

        !Mitja Added
        INTEGER howfound(2) ! (1) num found, (2) num not found
        INTEGER inthowfound(3) ! (1) num found, (2) num not found, (3) num multi found
        INTEGER lost
      END TYPE

      TYPE(CntType) :: cnt

      CONTAINS

!
! -----------------------------------------------------------------------------------------
!
      SUBROUTINE cntInit()
!
!     Set counters to zero
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE

      cnt%rTime=0.0D0
      cnt%iTime=0
      cnt%LastPID=-1
      cnt%cnp=0

      cnt%howfound = 0
      cnt%inthowfound = 0
      cnt%lost = 0

      END SUBROUTINE

!
! -----------------------------------------------------------------------------------------
!
      LOGICAL FUNCTION cntExport(i)
!
!     do i export
!
! -----------------------------------------------------------------------------------------
      INTEGER i

      cntExport=.FALSE.
      IF (i.GT.0) THEN
        IF (MODULO(cnt%iTime,i).EQ.0) cntExport=.TRUE.
      END IF

      END FUNCTION


      END MODULE
