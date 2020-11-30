!
!
!
      MODULE ParalelEnvironment
!
!     Include library headers
!  
      INCLUDE 'mpif.h'
!
! -----------------------------------------------------------------------------------------
!
      TYPE envType

        INTEGER          :: comm        ! communicator number
        INTEGER          :: nproc       ! number of processes
        INTEGER          :: mpr         ! my process number (zero based)
        INTEGER          :: myproc      ! my process number (one based)

      END TYPE

      TYPE(envType) :: env

      CONTAINS



! -----------------------------------------------------------------------------
      SUBROUTINE par_init()
!
!     $: Init parallel environment
!
! -----------------------------------------------------------------------------
      INTEGER ierr

      ierr=0

      env%comm = MPI_COMM_WORLD
      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(env%comm,env%nproc,ierr)
      CALL MPI_COMM_RANK(env%comm,env%mpr,ierr)
      env%myproc = env%mpr+1

      IF (ierr.NE.0) THEN
        Print *,"Could not init parallel environment!!"
        STOP
      END IF

      END SUBROUTINE

! -----------------------------------------------------------------------------
      INTEGER FUNCTION par_NextPID(LastPID)
!
!     $: Get ID of next processor
!
! -----------------------------------------------------------------------------      
      INTEGER LastPID

      par_NextPID = MOD(LastPID + 1,env%nproc)
      
      END FUNCTION


! -----------------------------------------------------------------------------
      LOGICAL FUNCTION par_AmIRankZero()
!
!     $: Check if rank 0
!
! -----------------------------------------------------------------------------
      
      IF (env%mpr.EQ.0) THEN
        par_AmIRankZero=.TRUE.
      ELSE
        par_AmIRankZero=.FALSE.
      END IF

      END FUNCTION

! -----------------------------------------------------------------------------
      LOGICAL FUNCTION par_IsMine(pid)
!
!     $: Check if rank 0
!
! -----------------------------------------------------------------------------
      INTEGER pid
      IF (env%mpr.EQ.pid) THEN
        par_IsMine=.TRUE.
      ELSE
        par_IsMine=.FALSE.
      END IF

      END FUNCTION


      END MODULE