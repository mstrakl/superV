!
! -----------------------------------------------------------------------------------------
!
      MODULE sorting
!
! -----------------------------------------------------------------------------------------
!
      implicit none
      public :: QsortC
      public :: bubles
      private :: Partition
    
      contains



! -----------------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE  QsortC(A, index)
!
! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing
!      
! -----------------------------------------------------------------------------------------
      REAL(8), intent(in out), dimension(:) :: A
      INTEGER, intent(in out), dimension(:) :: index
      integer :: iq
      
      IF (size(A) > 1) THEN
        CALL Partition(A, index, iq)
        CALL QsortC(A(:iq-1),index(:iq-1))
        CALL QsortC(A(iq:),index(iq:))
      END IF

      END SUBROUTINE

! -----------------------------------------------------------------------------------------    
      SUBROUTINE Partition(A, index, marker)
!      
! -----------------------------------------------------------------------------------------      
      INTEGER :: n
      real(8), intent(in out), dimension(:) :: A
      INTEGER, intent(in out), dimension(:) :: index
      integer, intent(out) :: marker
      integer :: i, j, itmp
      real(8) :: temp
      real(8) :: x      ! pivot point
      x = A(1)
      i= 0
      j= size(A) + 1
    
      do
         j = j-1
         do
            if (A(j) <= x) exit
            j = j-1
         end do
         i = i+1
         do
            if (A(i) >= x) exit
            i = i+1
         end do
         if (i < j) then
            ! exchange A(i) and A(j)
            temp = A(i)
            A(i) = A(j)
            A(j) = temp
            itmp = index(i)
            index(i)=index(j)
            index(j)=itmp
         elseif (i == j) then
            marker = i+1
            return
         else
            marker = i
            return
         endif
      end do
    
      END SUBROUTINE

!
!----------------------------------------------------------------------c
!                                                                      c
      subroutine bubles(nx,nix,lenix,x,ix,str)
! description:
! bublesort vector ix depending on x
! str= ... asc ... from small to big value
! str= ... des ... from big to small value
! call bubles(m%nbe,m%nbe,n,xy,ixy,"asc")
!----------------------------------------------------------------------c
!
         integer nx,nix,lenix,ix(nix)
         real(8) x(nx)
         REAL(8), ALLOCATABLE :: tmp(:)
         character(len=3):: str
! local
         integer i,j,k
!
!*** loop over position on ix
!
         if (str.eq.'asc') then
           do j=1,lenix-1
!
!*** loop over pairs
!
             do i=1,lenix-j
               if (x(ix(i)).gt.x(ix(i+1))) then
! change index
                 k=ix(i+1)
                 ix(i+1)=ix(i)
                 ix(i)=k
               endif
             enddo !i=j,lenix-1
           enddo !j=1,lenix-1
         else if  (str.eq.'des') then
           do j=1,lenix-1
!
!*** loop over pairs
!
             do i=1,lenix-j
               if (x(ix(i)).lt.x(ix(i+1))) then
! change index
                 k=ix(i+1)
                 ix(i+1)=ix(i)
                 ix(i)=k
               endif
             enddo !i=j,lenix-1
           enddo !j=1,lenix-1
         else
           stop 'bubles.f: wrong str'
         endif

!        sort the original vector
         ALLOCATE (tmp(nx))
         DO i=1,nx
            tmp(i)=x(ix(i))
         END DO
         DO i=1,nx
            x(i)=tmp(i)
         END DO
         DEALLOCATE(tmp)

!
         return
         end SUBROUTINE
   


      END MODULE
