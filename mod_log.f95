!
!
!
      MODULE logFile
!
! -----------------------------------------------------------------------------------------
!
      TYPE LogType
        INTEGER lun
        CHARACTER(100) FileName,IDname,IDversion,IDdate
        CHARACTER(39) StartTime,EndTime
      END TYPE

      TYPE(LogType) :: lg

      CONTAINS
! -----------------------------------------------------------------------------------------
      SUBROUTINE logInit()
!
!     Init log file
!
! -----------------------------------------------------------------------------------------
      USE ParalelEnvironment
      IMPLICIT NONE
      
      IF (par_AmIRankZero()) THEN
        lg%lun = 11
        WRITE (lg%FileName,'(A,A)') TRIM(lg%IDname),".log"
        CALL logGetDateTime(lg%StartTime)

        OPEN (lg%lun,FILE=TRIM(lg%FileName),STATUS='UNKNOWN',ERR = 10)

        CALL logWrite(TRIM(lg%IDname)//", v"//TRIM(lg%IDversion)//", "//TRIM(lg%IDdate))
        CALL logWrite(lg%StartTime)
        CALL logIntWrite("Number of processors:",env%nproc)
      END IF

      RETURN

10    CONTINUE
      WRITE (*,*) "Cant open log file!"
      STOP      

      END SUBROUTINE

! -----------------------------------------------------------------------------------------
      SUBROUTINE logWrite(text)
!
!     Write to log file
!
! -----------------------------------------------------------------------------------------
      USE ParalelEnvironment
      IMPLICIT NONE
      CHARACTER*(*) text
      CHARACTER(6) time      

      IF (par_AmIRankZero()) THEN
        CALL logGetTime(time)
        WRITE(lg%lun,'(A6,A)') time,text
        FLUSH(lg%lun)
      END IF

      END SUBROUTINE

! -----------------------------------------------------------------------------------------
      SUBROUTINE logIntWrite(text,I)
!
!     Write to log file
!
! -----------------------------------------------------------------------------------------
      USE ParalelEnvironment
      IMPLICIT NONE
      CHARACTER*(*) text
      CHARACTER(6) time
      INTEGER I

      IF (par_AmIRankZero()) THEN
        CALL logGetTime(time)
        WRITE(lg%lun,'(A6,A,I0)') time,text,I
        FLUSH(lg%lun)
      END IF

      END SUBROUTINE

! -----------------------------------------------------------------------------------------
      SUBROUTINE logPercent(i,n,p)
!
!     Write to log file
!
! -----------------------------------------------------------------------------------------
      USE ParalelEnvironment
      IMPLICIT NONE
      INTEGER i,n,v,p

      IF (par_AmIRankZero()) THEN
        v=MAX(n/(100/p),1)
        IF (MOD(i,v).EQ.0) THEN
          CALL logRealWrite("% ",100.0*i/DBLE(n),'(A6,A,F5.1)')
        END IF
      END IF

      END SUBROUTINE

! -----------------------------------------------------------------------------------------
      SUBROUTINE logRealWrite(text,R,fmt)
!
!     Write to log file
!
! -----------------------------------------------------------------------------------------
      USE ParalelEnvironment
      IMPLICIT NONE
      CHARACTER*(*) text,fmt
      CHARACTER(6) time      
      REAL(8) R

      IF (par_AmIRankZero()) THEN
        CALL logGetTime(time)
        WRITE(lg%lun,TRIM(fmt)) time,TRIM(text),R
        FLUSH(lg%lun)
      END IF

      END SUBROUTINE

! -----------------------------------------------------------------------------------------
      SUBROUTINE logVecWrite(text,R,n,fmt)
!
!     Write to log file
!
! -----------------------------------------------------------------------------------------
      USE ParalelEnvironment
      IMPLICIT NONE
      CHARACTER*(*) text,fmt
      CHARACTER(6) time   
      INTEGER i,n   
      REAL(8) R(n)

      IF (par_AmIRankZero()) THEN
        CALL logGetTime(time)
        WRITE(lg%lun,TRIM(fmt)) time,text,(R(i),i=1,n)
        FLUSH(lg%lun)
      END IF

      END SUBROUTINE



! -----------------------------------------------------------------------------------------
      SUBROUTINE logClose()
!
!     Close log file
!
! -----------------------------------------------------------------------------------------
      USE ParalelEnvironment
      IMPLICIT NONE

      IF (par_AmIRankZero()) THEN
        CALL logGetDateTime(lg%EndTime)
        CALL logWrite(lg%EndTime)

        CLOSE(lg%lun)
      END IF

      END SUBROUTINE

! -----------------------------------------------------------------------------
      SUBROUTINE logGetDateTime(cas)
!
!     $: writes date and time to cas
!
! -----------------------------------------------------------------------------            
      CHARACTER*50 D,T,Z
      INTEGER Value(8)
      CHARACTER*39 cas     
      
      CALL DATE_AND_TIME(D,T,Z,Value)
      WRITE (cas,33)  &
     &    'Date and time = ',Value(3),'.',Value(2),'.',Value(1),' ', &
         Value(5),':',Value(6),':',Value(7),'.',Value(8)
33    FORMAT (A16,I2,A1,I2,A1,I4,A1,I2,A1,I2,A1,I2,A1,I3)

      END SUBROUTINE   


! -----------------------------------------------------------------------------
      SUBROUTINE logGetTime(cas)
!
!     $: writes date and time to cas
!
! -----------------------------------------------------------------------------            
      CHARACTER*50 D,T,Z
      INTEGER Value(8)
      CHARACTER*6 cas     
      
      CALL DATE_AND_TIME(D,T,Z,Value)
      WRITE (cas,'(I2,A1,I2,1X)') Value(5),':',Value(6) 

      END SUBROUTINE    


      END MODULE
