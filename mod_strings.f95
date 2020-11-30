      MODULE String_Utility 
      IMPLICIT NONE 
      PRIVATE 
      PUBLIC :: StrLowCase 
      PUBLIC :: isNumeric

   ! List of character for case conversion 
      CHARACTER(*), PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz' 
      CHARACTER(*), PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 

      CONTAINS 

      FUNCTION StrLowCase( Input_String ) RESULT( Output_String ) 
        CHARACTER(*), INTENT(IN)     :: Input_String 
        CHARACTER(LEN(Input_String)) :: Output_String 
        INTEGER :: i, n 

        ! Copy input string 
        Output_String = Input_String 

        ! Convert case character by character 
        DO i = 1, LEN(Output_String) 
          n = INDEX(UPPER_CASE, Output_String(i:i)) 
          IF ( n /= 0 ) Output_String(i:i) = LOWER_CASE(n:n) 
        END DO 
      END FUNCTION StrLowCase 
!______________________________________________________________________!
!______________________________________________________________________!

      FUNCTION isNumeric(string)
        IMPLICIT NONE
        CHARACTER(len=*), INTENT(IN) :: string
        LOGICAL :: isNumeric
        REAL :: x
        INTEGER :: e,n
        CHARACTER(len=12) :: fmt
        
        n=LEN_TRIM(string)
        WRITE(fmt,'("(F",I0,".0)")') n
        READ(string,fmt,IOSTAT=e) x
        isNumeric = e == 0
      END FUNCTION isNumeric 
      
      !OneLine=ADJUSTL(OneLine)  ! remove leadin spaces
      !READ(OneLine,*) tmp(i) ! read the first number in the line
      !ipos = SCAN(TRIM(OneLine)," ") ! find next space

      END MODULE String_Utility 