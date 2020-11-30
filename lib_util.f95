      subroutine create_directory( newDirPath )
      ! Author:  Jess Vriesema
      ! Date:    Spring 2011
      ! Purpose: Creates a directory at ./newDirPath

      implicit none

      character(len=*), intent(in) :: newDirPath
      character(len=256)           :: mkdirCmd
      logical                      :: dirExists

      ! Check if the directory exists first
      inquire( file=trim(newDirPath)//'/.', exist=dirExists )  ! Works with gfortran, but not ifort
      !inquire( directory=newDirPath, exist=dirExists )         ! Works with ifort, but not gfortran


      if (dirExists) then
!      write (*,*) "Directory already exists: '"//trim(newDirPath)//"'"
      else
        mkdirCmd = 'mkdir -p '//trim(newDirPath)
!        write(*,'(a)') "Creating new directory: '"//trim(mkdirCmd)//"'"
        call system( mkdirCmd )
      endif
      end subroutine create_directory

      
! *********************************************************************
! **                                                                 **
! ** Squeeze multiple blanks into single blank                       **
! **                                                                 **
! *********************************************************************
!---------------------------------------------------------------------C
      SUBROUTINE sqblnk(lun,line)
!---------------------------------------------------------------------C
      CHARACTER*(*) line
      LOGICAL flag
      INTEGER i,ii,j,lnblnk,lun

      j=1
      flag=.FALSE.
      DO i=1,LNBLNK(line)+1
        IF (line(i:i).NE.' ' .AND. .NOT.flag) THEN
          flag=.TRUE.
          ii=i
        ELSE IF (line(i:i).EQ.' ' .AND. flag) THEN
          flag=.FALSE.
          line(j:j+i-ii)=line(ii:i)
          j=j+i-ii+1
        END IF
      END DO
      WRITE(lun,'(A)') line(1:j-2)
      RETURN
      END
      INTEGER FUNCTION LNBLNK (string)
!
!     LNBLNK returns the index of the last non-blank character in string
!
      CHARACTER*(*) string
      INTEGER i

      DO i=LEN(string),1,-1
        IF (string(i:i).NE.' ') THEN
          lnblnk=i
          RETURN
        END IF
      END DO
      lnblnk=0
      RETURN
      END           

!______________________________________________________________________C
!______________________________________________________________________C
      SUBROUTINE rOneTL(lun,OneLine)
!     _    ___ _    _ 
!     Read One Text Line 
!
!______________________________________________________________________C
!     Returns the first nonempty text line in file LUN, which does not
!     include the # character. If end of file is encoutered, it returns EOF
      CHARACTER*(*) OneLine
      INTEGER lun,i
  
10    READ(lun,'(A)',END=20) OneLine  

!     Check if line is empty
      IF (len_trim(OneLine).EQ.0) GOTO 10

!     Check if line contains # character
      DO i=1,len_trim(OneLine)
        IF (OneLine(i:i).EQ.'#') GOTO 10
      ENDDO

      RETURN

20    OneLine='EOF'
      END           
   
!______________________________________________________________________!
!______________________________________________________________________!
      SUBROUTINE StopProgram(ierr)
!     _    ___ _    _ 
!     Read One Text Line 
!
      USE logFile
      USE ParalelEnvironment
      IMPLICIT NONE 

      INTEGER ierr

      IF (ierr.NE.0) THEN
        CALL logWrite("Finished with a error!")
      ELSE
        CALL logWrite("Finished sucessfully!")
      END IF

      CALL logClose()
      CALL MPI_FINALIZE(env%comm)

      STOP

      END


!----------------------------------------------------------------------C
!----------------------------------------------------------------------c
!                                                                      c
      SUBROUTINE rdvec(ifr,nnx,vec)
!                                                                      c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------C
!.......................................................................
!..                                                                   ..
!..   REad VECtor using MAXSIZE chunks                                ..
!..   --   ---                                                        ..
!.......................................................................
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      REAL*8  vec(nnx)
      PARAMETER (maxsize=8192)

      nblock=INT(nnx/maxsize)
      DO j=1,nblock
        k=maxsize*(j-1)
        READ(ifr) (vec(i),i=k+1,k+maxsize)
      END DO
      k=maxsize*nblock
      READ(ifr) (vec(i),i=k+1,nnx)
      RETURN
      END
!----------------------------------------------------------------------C
!----------------------------------------------------------------------c
!                                                                      c
      SUBROUTINE wrvec(ifr,nnx,vec)
!                                                                      c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------C
!.......................................................................
!..                                                                   ..
!..   WRite VECtor using MAXSIZE chunks                               ..
!..   --    ---                                                       ..
!.......................................................................
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      REAL*8  vec(nnx)
      PARAMETER (maxsize=8192)
!....
      nblock=INT(nnx/maxsize)
      DO j=1,nblock
        k=maxsize*(j-1)
        WRITE(ifr) (vec(i),i=k+1,k+maxsize)
      END DO
      k=maxsize*nblock
      WRITE(ifr) (vec(i),i=k+1,nnx)
      RETURN
      END

!______________________________________________________________________C
!______________________________________________________________________C
      SUBROUTINE WrMat(mat,nrow,ncol,io)
!        __    _      ___
!        Write Single Matrix
!______________________________________________________________________C
!______________________________________________________________________C
      INTEGER nrow,ncol,io,j
      REAL(8) mat(nrow,ncol)

      DO j=1,ncol
        CALL wrvec(io,nrow,mat(1,j))
      END DO

      END
!______________________________________________________________________C
!______________________________________________________________________C
      SUBROUTINE RdMat(mat,nrow,ncol,io)
!        _  _ _      ___
!        Read Single Matrix
!______________________________________________________________________C
!______________________________________________________________________C
      INTEGER nrow,ncol,io,j
      REAL(8) mat(nrow,ncol)

      DO j=1,ncol
        CALL rdvec(io,nrow,mat(1,j))
      END DO

      END


!______________________________________________________________________C
!______________________________________________________________________C
      SUBROUTINE WrIMat(mat,nrow,ncol,io)
!        __    _      ___
!        Write Single Matrix
!______________________________________________________________________C
!______________________________________________________________________C
      INTEGER nrow,ncol,io,j
      INTEGER mat(nrow,ncol)

      DO j=1,ncol
        CALL wrIvec(io,nrow,mat(1,j))
      END DO

      END
!______________________________________________________________________C
!______________________________________________________________________C
      SUBROUTINE RdIMat(mat,nrow,ncol,io)
!        _  _ _      ___
!        Read Single Matrix
!______________________________________________________________________C
!______________________________________________________________________C
      INTEGER nrow,ncol,io,j
      INTEGER mat(nrow,ncol)

      DO j=1,ncol
        CALL rdIvec(io,nrow,mat(1,j))
      END DO

      END

!----------------------------------------------------------------------C
!----------------------------------------------------------------------c
!                                                                      c
      SUBROUTINE rdIvec(ifr,nnx,vec)
!                                                                      c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------C
!.......................................................................
!..                                                                   ..
!..   REad VECtor using MAXSIZE chunks                                ..
!..   --   ---                                                        ..
!.......................................................................
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      INTEGER  vec(nnx)
      PARAMETER (maxsize=8192)

      nblock=INT(nnx/maxsize)
      DO j=1,nblock
        k=maxsize*(j-1)
        READ(ifr) (vec(i),i=k+1,k+maxsize)
      END DO
      k=maxsize*nblock
      READ(ifr) (vec(i),i=k+1,nnx)
      RETURN
      END
!----------------------------------------------------------------------C
!----------------------------------------------------------------------c
!                                                                      c
      SUBROUTINE wrIvec(ifr,nnx,vec)
!                                                                      c
!----------------------------------------------------------------------c
!----------------------------------------------------------------------C
!.......................................................................
!..                                                                   ..
!..   WRite VECtor using MAXSIZE chunks                               ..
!..   --    ---                                                       ..
!.......................................................................
      INTEGER ifr,nnx,maxsize,nblock,i,j,k
      INTEGER  vec(nnx)
      PARAMETER (maxsize=8192)
!....
      nblock=INT(nnx/maxsize)
      DO j=1,nblock
        k=maxsize*(j-1)
        WRITE(ifr) (vec(i),i=k+1,k+maxsize)
      END DO
      k=maxsize*nblock
      WRITE(ifr) (vec(i),i=k+1,nnx)
      RETURN
      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE CalPolygonCentroid(x,npnts,xc)
!
!     Calculate polygon centroid by geometric decomposition onto triangles
!
!          e ---d 
!         /      \
!        /       c
!       /       /
!      a -----b
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER npnts,i
      REAL(8) x(npnts,3),xc(3)
      !REAL(8) a(3),b(3),n(3)
      !REAL(8) nn(npnts,3)

      xc = 0.0D0

      DO i = 1,npnts

        xc = xc + x(i,:)

      END DO

      IF (npnts.GT.0) THEN
        xc = (1.0D0/npnts) * xc
      ELSE
        xc = 0.0D0
        PRINT *, "ERROR :: CalPolygonCentroid :: npnts = 0 !!!"
      END IF
   
      END      

! -----------------------------------------------------------------------------------------
      SUBROUTINE CalPolygonNormal(x,npnts,n)
!
!     Calculate normal from multiple poins with Newell's method
!
!          e ---d 
!         /      \
!        /       c
!       /       /
!      a -----b
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER npnts,i
      REAL(8) x(npnts,3),xc(3)
      REAL(8) a(3),b(3),n(3),vp(3)

      n = 0.0D0

      IF (npnts.GT.3) THEN

        DO i = 1,npnts

          IF (i.EQ.1) THEN

            a = x(i+1,:) - x(i,:)
            b = x(npnts,:) - x(i,:)

          ELSE IF (i.EQ.npnts) THEN

            a = x(1,:) - x(i,:)
            b = x(i-1,:) - x(i,:)

          ELSE

            a = x(i+1,:) - x(i,:)
            b = x(i-1,:) - x(i,:)

          END IF

          CALL VectorProduct(a,b,vp)
          n = n + vp

        END DO

        CALL Normalize(n)

      ELSE

      CALL CalNormal3p( x(1,:),x(2,:),x(3,:),n )

      END IF


   
      END   

! -----------------------------------------------------------------------------------------
      SUBROUTINE CalNormal3p(a,b,c,n)
!
!     Calculate normal from three poins
!
!          c
!         /
!        /
!      a ---- b
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) a(3),b(3),c(3),ba(3),ca(3),n(3)

      ba=b-a
      ca=c-a

      CALL VectorProduct(ba,ca,n)
      CALL Normalize(n)

      END

      
! -----------------------------------------------------------------------------------------
      SUBROUTINE CalNormalMultiPoint(pnts,npnts,center,normal)
!
!     Calculate normal from multiple poins
!
!          e ---d 
!         /      \
!        /       c
!       /       /
!      a -----b
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER npnts,i
      REAL(8) pnts(npnts,3),center(3),normal(3)
      REAL(8) a(3),b(3),n(3)
      REAL(8) nn(npnts,3)

      DO i = 1,(npnts-1)

            a = pnts(i,:)-center(:)
            b = pnts(i+1,:)-center(:)

            CALL VectorProduct(a,b,n)
            CALL Normalize(n)

            nn(i,:) = n
      END DO

      i=npnts

      a = pnts(i,:)-center(:)
      b = pnts(1,:)-center(:)

      CALL VectorProduct(a,b,n)
      CALL Normalize(n)

      nn(i,:) = n

      normal(1) = 1.0D0*sum(nn(:,1))/npnts
      normal(2) = 1.0D0*sum(nn(:,2))/npnts
      normal(3) = 1.0D0*sum(nn(:,3))/npnts
      CALL Normalize(normal)

      END
! -----------------------------------------------------------------------------------------
      SUBROUTINE vecLen(v,d)
!
!     Vector length
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) v(3),d

      d = SQRT(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE dist2P (a,b,dsq)
!
!     Distance between 2 points
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) a(3),b(3),dsq

      dsq = SQRT( ( a(1)-b(1) )**2 + ( a(2)-b(2) )**2 + ( a(3)-b(3) )**2 )

      END
      
! -----------------------------------------------------------------------------------------
      SUBROUTINE DotProduct(a,b,r)
!
!     Calculate dot product
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) a(3),b(3),r

      r = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE VectorProduct(a,b,r)
!
!     Calculate vector product
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) a(3),b(3),r(3)

      r(1) = + a(2)*b(3) - a(3)*b(2)
      r(2) = - a(1)*b(3) + a(3)*b(1)
      r(3) = + a(1)*b(2) - a(2)*b(1)

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE Normalize(n)
!
!     Normalize vector
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) n(3),d


      d = SQRT(n(1)*n(1)+n(2)*n(2)+n(3)*n(3))   

      IF (d.GT.0.0D0) THEN
        n = n/d
      ELSE
        n = n
        PRINT *, "ERROR :: Normalize :: zero-length !!!"
        !CALL stopProgram(1)
      END IF

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE CalNormal2vec(a,b,n)
!
!     Calculate normal from 2 vectors
!
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) a(3),b(3),n(3)
      
      CALL VectorProduct(a,b,n)
      CALL Normalize(n)

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE vecAngle(a,b,phi)
!
!     Angle between vectors from dot product
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) a(3),b(3),phi
      REAL(8) dp,alen,blen

      CALL DotProduct(a,b,dp)
      CALL vecLen(a,alen)
      CALL vecLen(b,blen)

      phi = ACOS(dp/(alen*blen))

      END   
      

! -----------------------------------------------------------------------------------------
      SUBROUTINE GetBarycentricInTri(p,x,c,n,fb,ierr)
!
!     
! -----------------------------------------------------------------------------------------
      
      IMPLICIT NONE
      REAL(8) p(3),x(3,3),porg(3)
      REAL(8) fb(3)
      INTEGER ierr

      REAL(8) vp(3),Abary(3),A,small

      REAL(8) c(3),n(3),dist,dist2

      small = 1.0E-5
      ierr = 0

      porg = p

      !Perpendicular point - plane distance
      CALL DotProduct( p-c, n, dist )

      !Project point to a plane
      p = p - dist*n

      !Perpendicular point - plane distance
      CALL DotProduct( p-c, n, dist2 )

      !Area x1 - x2 - p
      CALL VectorProduct( x(1,:)-p, x(2,:)-p, vp )
      CALL vecLen(vp,Abary(3))

      !Area x2 - x3 - p
      CALL VectorProduct( x(2,:)-p, x(3,:)-p, vp )
      CALL vecLen(vp,Abary(1))

      !Area x3 - x1 - p
      CALL VectorProduct( x(3,:)-p, x(1,:)-p, vp )
      CALL vecLen(vp,Abary(2))

      !Area triangle
      CALL VectorProduct( x(2,:)-x(1,:), x(3,:)-x(1,:), vp )
      CALL vecLen(vp,A)

      fb(1) = Abary(1) / A
      fb(2) = Abary(2) / A
      fb(3) = Abary(3) / A

      IF (fb(1)+fb(2)+fb(3).GT.1.0D0+small) ierr = 1
      IF (fb(1)+fb(2)+fb(3).LT.1.0D0-small) ierr = 1

      IF (dist2.GT.1.0E-10 .AND. ierr.EQ.0) THEN

            !WRITE(*,'(A4,F12.6)') "A1=",Abary(1)
            !WRITE(*,'(A4,F12.6)') "A2=",Abary(2)
            !WRITE(*,'(A4,F12.6)') "A3=",Abary(3)
            !WRITE(*,'(A4,F12.6)') "A=",A

            WRITE(*,'(A4,3F20.12)') "x1=",x(1,:)
            WRITE(*,'(A4,3F20.12)') "x2=",x(2,:)
            WRITE(*,'(A4,3F20.12)') "x3=",x(3,:)
            WRITE(*,'(A4,3F20.12)') "fc=",c
            WRITE(*,'(A4,3F20.12)') "fn=",n
            WRITE(*,'(A5,3F20.12)') "porg=",porg
            WRITE(*,'(A6,3F20.12)') "pproj=",p
            WRITE(*,'(A4,3F20.12,G20.12)') "fbs=",fb,sum(fb)

            WRITE(*,'(A5,G20.12)') "dist=",dist
            WRITE(*,'(A6,G20.12)') "dist2=",dist2   
            
            CALL BreakPoint()

      END IF




!      IF(ierr.EQ.0)  THEN
!      
!            WRITE(*,'(A4,F12.6)') "A1=",Abary(1)
!            WRITE(*,'(A4,F12.6)') "A2=",Abary(2)
!            WRITE(*,'(A4,F12.6)') "A3=",Abary(3)
!            WRITE(*,'(A4,F12.6)') "A=",A
!
!            WRITE(*,'(A4,3F12.6)') "x1=",x(1,:)
!            WRITE(*,'(A4,3F12.6)') "x2=",x(2,:)
!            WRITE(*,'(A4,3F12.6)') "x3=",x(3,:)
!            WRITE(*,'(A4,3F12.6)') "fc=",c
!            WRITE(*,'(A4,3F12.6)') "fn=",n
!
!            WRITE(*,'(A4,3F12.6)') "p=",p
!
!            WRITE(*,'(A4,G14.6)') "dp=",dist
!            WRITE(*,'(A4,3F12.6,G14.6)') "fbs=",fb,sum(fb)
!
!
!            !IF (abs(dist).GT.1E-6) THEN
!              !CALL WriteVTKPoint('point.vtk',p)
!              !CALL WriteVTKTriangle('triangle.vtk',x(1,:),x(2,:),x(3,:),c,n)
!            !END IF
!      END IF


      END  

 ! -----------------------------------------------------------------------------------------
      SUBROUTINE GetBarycentricInTriDebug(p,x,c,n,fb,ierr)
!
!     
! -----------------------------------------------------------------------------------------
      
      IMPLICIT NONE
      REAL(8) p(3),x(3,3),porg(3)
      REAL(8) fb(3)
      INTEGER ierr

      REAL(8) vp(3),Abary(3),A,small

      REAL(8) c(3),n(3),dist,dist2

      small = 1.0E-8
      ierr = 0

      porg = p

      !Perpendicular point - plane distance
      CALL DotProduct( p-c, n, dist )

      !Project point to a plane
      p = p - dist*n

      !Perpendicular point - plane distance
      CALL DotProduct( p-c, n, dist2 )

      !Area x1 - x2 - p
      CALL VectorProduct( x(1,:)-p, x(2,:)-p, vp )
      CALL vecLen(vp,Abary(3))

      !Area x2 - x3 - p
      CALL VectorProduct( x(2,:)-p, x(3,:)-p, vp )
      CALL vecLen(vp,Abary(1))

      !Area x3 - x1 - p
      CALL VectorProduct( x(3,:)-p, x(1,:)-p, vp )
      CALL vecLen(vp,Abary(2))

      !Area triangle
      CALL VectorProduct( x(2,:)-x(1,:), x(3,:)-x(1,:), vp )
      CALL vecLen(vp,A)

      fb(1) = Abary(1) / A
      fb(2) = Abary(2) / A
      fb(3) = Abary(3) / A

      IF (fb(1)+fb(2)+fb(3).GT.1.0D0+small) ierr = 1
      IF (fb(1)+fb(2)+fb(3).LT.1.0D0-small) ierr = 1

      IF (dist.GT.1.0E-10 .AND. ierr.EQ.0) THEN

            !WRITE(*,'(A4,F12.6)') "A1=",Abary(1)
            !WRITE(*,'(A4,F12.6)') "A2=",Abary(2)
            !WRITE(*,'(A4,F12.6)') "A3=",Abary(3)
            !WRITE(*,'(A4,F12.6)') "A=",A

            WRITE(*,'(A4,3F20.12)') "x1=",x(1,:)
            WRITE(*,'(A4,3F20.12)') "x2=",x(2,:)
            WRITE(*,'(A4,3F20.12)') "x3=",x(3,:)
            WRITE(*,'(A4,3F20.12)') "fc=",c
            WRITE(*,'(A4,3F20.12)') "fn=",n
            WRITE(*,'(A5,3F20.12)') "porg=",porg
            WRITE(*,'(A6,3F20.12)') "pproj=",p
            WRITE(*,'(A4,3F20.12,G20.12)') "fbs=",fb,sum(fb)

            WRITE(*,'(A5,G20.12)') "dist=",dist
            WRITE(*,'(A6,G20.12)') "dist2=",dist2   
            
            CALL BreakPoint()

      END IF




!      IF(ierr.EQ.0)  THEN
!      
!            WRITE(*,'(A4,F12.6)') "A1=",Abary(1)
!            WRITE(*,'(A4,F12.6)') "A2=",Abary(2)
!            WRITE(*,'(A4,F12.6)') "A3=",Abary(3)
!            WRITE(*,'(A4,F12.6)') "A=",A
!
!            WRITE(*,'(A4,3F12.6)') "x1=",x(1,:)
!            WRITE(*,'(A4,3F12.6)') "x2=",x(2,:)
!            WRITE(*,'(A4,3F12.6)') "x3=",x(3,:)
!            WRITE(*,'(A4,3F12.6)') "fc=",c
!            WRITE(*,'(A4,3F12.6)') "fn=",n
!
!            WRITE(*,'(A4,3F12.6)') "p=",p
!
!            WRITE(*,'(A4,G14.6)') "dp=",dist
!            WRITE(*,'(A4,3F12.6,G14.6)') "fbs=",fb,sum(fb)
!
!
!            !IF (abs(dist).GT.1E-6) THEN
!              !CALL WriteVTKPoint('point.vtk',p)
!              !CALL WriteVTKTriangle('triangle.vtk',x(1,:),x(2,:),x(3,:),c,n)
!            !END IF
!      END IF


      END      

! -----------------------------------------------------------------------------------------
      SUBROUTINE Breakpoint()
!
!     Compilation breakpoint
! -----------------------------------------------------------------------------------------
      
      use, intrinsic:: iso_fortran_env, only: stdin=>input_unit
      
      IMPLICIT NONE

      print *, 'Breakpoint, press Enter.'
      read(stdin,*)

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE makeTriangles(con,n,ibest,tricon)
!
! -----------------------------------------------------------------------------------------

      IMPLICIT NONE 
      !INTEGER n,nTri,nCon,i,j,k,l,ibest
      INTEGER i,j,ibest,n,nTri
      INTEGER con(n),newcon(n),tricon(n-2,3)
      !REAL(8) pnts(n,3)

      DO i = 1,n
        j=ibest+(i-1)
        IF (j .GT. n) j = j-n
          newcon(i) = con(j)	
          !WRITE(*,'(I4,3F10.5)') con(i),pnts(i,:)
      END DO

      nTri = n-2

      DO i = 1,nTri

        tricon(i,1) = newcon(1)
        tricon(i,2) = newcon(i+1)
        tricon(i,3) = newcon(i+2)

      END DO

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE calcAngles(pnts,n,ibest)
!	    Calculates which vertex forms best quality triangles from face
!
! -----------------------------------------------------------------------------------------

      IMPLICIT NONE 
      INTEGER n,i,j,k,l,ibest
      REAL(8) pnts(n,3),vecs(n-1,4),a(3),b(3),dp,alen,blen,cosA(n-2),minCos(n),maxMinCos

      ibest=1

      DO i = 1,n

        k = i
        DO j = 1,n-1

          k=k+1
          IF (k .GT. n) k = k-n
          vecs(j,1:3) = pnts(k,:)-pnts(i,:)
          vecs(j,4) = k

        END DO

        DO j = 1,n-2

          a = vecs(j,1:3)
          b = vecs(j+1,1:3)
          CALL DotProduct(a,b,dp)
          CALL vecLen(a,alen)
          CALL vecLen(b,blen)
          cosA(j)=dp/(alen*blen)

        END DO

        minCos(i) = minval(cosA)

      END DO	

      !WRITE(*,'(A10,6F10.5)') "minCos=",minCos

      ibest = MAXLOC(minCos,1)

      !maxMinCos = minCos(ibest)
      !PRINT *, ibest, maxMinCos
      !PRINT *, "GOING OUT"


      END




! -----------------------------------------------------------------------------------------
      LOGICAL FUNCTION ComaprePoints(a,b,eps)
!
!     Are point the same?
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) a(3),b(3),eps
      LOGICAL result

      ComaprePoints = .FALSE.

      IF ( (ABS(a(1)-b(1))<eps) .AND. (ABS(a(2)-b(2))<eps) .AND. (ABS(a(3)-b(3))<eps) ) THEN
            ComaprePoints = .TRUE.
      END IF

      END

! -----------------------------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION Distance2(a,b)
!
!     square of distance
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8) a(3),b(3)

      Distance2 = ( a(1) - b(1) )**2 + (a(2)-b(2) )**2 + (a(3)-b(3) )**2

      END

! -----------------------------------------------------------------------------------------
      INTEGER FUNCTION FindLoc(array,len,val)
!
!
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE

      INTEGER len,val
      INTEGER array(len)
      INTEGER i

      FindLoc = 0

      DO i = 1,len

            IF (array(i).EQ.val) THEN
                  FindLoc = i
                  RETURN
            END IF

      END DO    

      END
