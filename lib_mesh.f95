! -----------------------------------------------------------------------------------------
      SUBROUTINE GetVTKmesh(m,fname,fname2)
!
!     Read the mesh information from the VTK file
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile  
      USE ParalelEnvironment          
      IMPLICIT NONE 
      CHARACTER*(*) fname !VTK results path and filename
      CHARACTER*(*) fname2 !foamMesh path
      CHARACTER(255) bFname
      CHARACTER(255) OneLine
      TYPE(meshType) :: m  ! mesh data structure
      INTEGER ierr,hex,tet,pri,pyr,poly,i,j

      REAL t1,t2
!
!     Set element local coordinates
!
      CALL SetElementLC(m)

!
!     Read mesh points and get number of elements from VTK file
!
      CALL ReadMeshPointsVTK(m,fname)
      PRINT *, "Finished ReadMeshPointsVTK"
!
!     Read mesh files from OF "constant/polyMesh/ dir 
!     -used as supporting files for mesh elements and faces addresing
!        
      CALL ReadFoamMesh(m,fname2)
      PRINT *, "Finished ReadFoamMesh"      
!
!     Construct elements from foam faces, owners and neighbours 
! 
      CALL ConstructFoamElementsFromFaces(m)
      PRINT *, "Finished ConstructFoamElementsFromFaces"
!
!     Read mesh files from OF "constant/polyMesh/ dir 
!     -used as supporting files for mesh elements and faces addresing
!        
      CALL TriangulateFoamFaces(m)  
      PRINT *, "Finished TriangulateFoamFaces"
!
!     Check foam faces triangulations for possible zerovolume problem
!
      CALL CheckFoamTriangulation(m)
!
!     Construct elements starting from decomposed tri-faces
! 
      !CALL ConstructFoamElementsFromTriFaces(m)     
!
!     Calculate volumetric center of foam elements
      !CALL CalFoamElementsVolCenter(m,1)        
!
!     Check elements -> to baš ne funkcionira...
! 
      !CALL CheckFoamMeshElements(m)
!
!     Creat tetrahedrons from face triangles and element vol centers
! 
      CALL CreateFoamDecomposedTets(m)
      PRINT *, "Finished CreateFoamDecomposedTets"
!
!     Calculate face / element centers and normals 
!        
      CALL CalculateCentersNormals(m,1)
      PRINT *, "Finished CalculateCentersNormals"
!
!     Reached end of mesh reading
! 
      PRINT *, "Finished with mesh reading"
      !CALL stopProgram(1)
!
!     Remenber mesh name
!
      m%FullName=fname
!
!     Generate mesh report
!
      CALL logWrite    ("")
      CALL logWrite    ("Mesh stats")
      CALL logWrite    ("Mesh name (ASCII)  = "//TRIM(fname))
      CALL logWrite    ("Mesh name (BIN)    = "//TRIM(bFname))
      CALL logIntWrite ("Number of nodes    = ",m%nnodes)
      CALL logIntWrite ("Number of elements = ",m%nelem)

      CALL logWrite    ("Mesh extents")
      CALL logRealWrite("              xmin = ",m%xmin,'(A6,A,F10.5)')
      CALL logRealWrite("              xmax = ",m%xmax,'(A6,A,F10.5)')
      CALL logRealWrite("              ymin = ",m%ymin,'(A6,A,F10.5)')
      CALL logRealWrite("              ymax = ",m%ymax,'(A6,A,F10.5)')
      CALL logRealWrite("              zmin = ",m%zmin,'(A6,A,F10.5)')
      CALL logRealWrite("              zmax = ",m%zmax,'(A6,A,F10.5)')

      END
   
! -----------------------------------------------------------------------------------------
      SUBROUTINE ReadFoamMesh(m,fname)
!     
!     Read face connectivity from foam file 
! -----------------------------------------------------------------------------------------
      USE mesh
      USE String_Utility
      USE logFile   

      IMPLICIT NONE 
      CHARACTER KeyWord*64,OneLine*255,foamfile*50,dummy*64,PrevLine*64
      CHARACTER*(*) fname 

      TYPE(meshType) :: m
      INTEGER lun,i,j,k,nentries,ipos,br1,br2


      REAL t1,t2

      CALL CPU_TIME(t1)

! -----------------------------------------------------------------------------------------
!     Read mesh faces  
! -----------------------------------------------------------------------------------------      

      lun = 13
      foamfile = "/faces"
      OPEN (lun,FILE=TRIM(fname) // TRIM(foamfile),ERR=10,STATUS='OLD')

      i=0
      j=0
      m%nfac=0

      OneLine="NOT EOF"
      DO WHILE (OneLine(1:3).NE.'EOF')

        i=i+1
        CALL rOneTL(lun,OneLine)
        READ (OneLine,*) KeyWord

        !If num entries > 0 (successfull number read) 
        IF(m%nfac.GT.0) THEN

          !Check first and second bracket if they exist
          br1=SCAN(TRIM(OneLine),"(")
          br2=SCAN(TRIM(OneLine),")")

          !Procees only if both brackets present -> valid line
          IF (br1.GT.0.AND.br2.GT.0) THEN

            j=j+1

            !Read vector size -> this represents number of face vertices
            READ( OneLine(1:br1-1),* ) m%f(j)%nVertex
            !PRINT *, "vs=",m%f(j)%nVertex

            !Trim away size int + leading and trailind bracket
            OneLine=OneLine(br1+1:br2-1)
            !PRINT *, "ol=",OneLine

            ALLOCATE( m%f(j)%con( m%f(j)%nVertex ) )
            m%f(j)%con(:) = -1

            DO k = 1,m%f(j)%nVertex

              !Get vec element
              READ (OneLine,*) m%f(j)%con(k)

              !add +1 for difference in count start (Fortran starts with 1, OpenFOAM with 0)
              m%f(j)%con(k) = m%f(j)%con(k) + 1

              !PRINT *, "vec(",k,")=",m%f(j)%con(k)

              !Account for spaces and trim line
              ipos = SCAN(TRIM(OneLine)," ") ! find next space
              OneLine=OneLine(ipos+1:LEN_TRIM(OneLine))

            END DO

            !PRINT *, "m%f(",j,")%con=",m%f(j)%con(:)

          END IF

        END IF

        !Read first number in file -> defines num entries
        IF (isNumeric(KeyWord).AND.m%nfac.EQ.0) THEN
          READ(KeyWord,*) m%nfac

          !Allocate space for number of faces
          ALLOCATE(m%f(m%nfac))
        END IF

      END DO

      IF (j.EQ.m%nfac) THEN
        PRINT *, "Foam face read OK; nfaces=",j
      ELSE
        PRINT *, "Foam face read m%nfac error; nfaces=",j
      END IF

      CLOSE(lun)

      CALL CPU_TIME(t2)
      PRINT *, "ReadFoamMesh elapsed time :: Faces:",t2-t1

! -----------------------------------------------------------------------------------------
!     Read mesh boundary  
! -----------------------------------------------------------------------------------------      

      lun = 13
      foamfile = "/boundary"
      OPEN (lun,FILE=TRIM(fname) // TRIM(foamfile),ERR=10,STATUS='OLD')

      i=0
      j=0
      PrevLine = "FALSE"

      OneLine="NOT EOF"
      DO WHILE (OneLine(1:3).NE.'EOF')

        i=i+1
        CALL rOneTL(lun,OneLine)
        READ (OneLine,*) KeyWord

        !If num entries > 0 (successfull number read) 
        IF(m%nbound.GT.0) THEN


          IF ( SCAN(PrevLine,"(").GT.0 .OR. SCAN(PrevLine,"}").GT.0 ) THEN

            !Skip if non-letters
            IF ( SCAN(OneLine,")").EQ.0 .AND. SCAN(OneLine,"}").EQ.0 ) THEN

              j = j + 1
              m%bound(j)%name = TRIM(KeyWord)

              DO WHILE ( SCAN(OneLine,"}").EQ.0 )

                CALL rOneTL(lun,OneLine)
                READ (OneLine,*) KeyWord

                IF (StrLowCase(TRIM(KeyWord)).EQ.'type') THEN
                  READ(OneLine,*) dummy,KeyWord
                    m%bound(j)%type = 0 !default
                  IF (TRIM(KeyWord).EQ.'wall') THEN
                    m%bound(j)%type = 1
                  ELSE IF (TRIM(KeyWord).EQ.'patch') THEN
                    m%bound(j)%type = 2
                  END IF

                ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'nfaces') THEN
                  READ(OneLine,*) dummy,m%bound(j)%nFaces
                ELSE IF (StrLowCase(TRIM(KeyWord)).EQ.'startface') THEN
                  READ(OneLine,*) dummy,m%bound(j)%startFace 
                END IF

              END DO

              !PRINT *, "Name(",j,")=",m%bound(j)%name
              !PRINT *, "type=",m%bound(j)%type
              !PRINT *, "nf=",m%bound(j)%nFaces
              !PRINT *, "sf=",m%bound(j)%startFace 
              !CALL BreakPoint()

            END IF

          END IF

          IF (LEN(TRIM(OneLine)).GT.0) PrevLine = TRIM(OneLine)

        END IF

        !Read first number in file -> defines num entries
        IF (isNumeric(KeyWord).AND.m%nbound.EQ.0) THEN

          !Allocate space for boundaries
          READ(KeyWord,*) m%nbound
          ALLOCATE(m%bound(m%nbound))

        END IF

      END DO

      CLOSE(lun)

      CALL CPU_TIME(t2)
      PRINT *, "ReadFoamMesh elapsed time :: Boundary:",t2-t1


! -----------------------------------------------------------------------------------------
!     Read face owners 
! -----------------------------------------------------------------------------------------      

      lun = 13
      foamfile = "/owner"
      OPEN (lun,FILE=TRIM(fname) // TRIM(foamfile),ERR=10,STATUS='OLD')

      i=0
      j=0
      nentries=0

      OneLine="NOT EOF"
      DO WHILE (OneLine(1:3).NE.'EOF')

        i=i+1
        CALL rOneTL(lun,OneLine)
        READ (OneLine,*) KeyWord

        !PRINT *, "ol=",TRIM(OneLine)

        !If num entries > 0 (successfull number read) 
        IF(nentries.GT.0) THEN

          !Check if leading or terminating bracket is found
          br1=SCAN(TRIM(OneLine),"(")
          br2=SCAN(TRIM(OneLine),")")

          !When br2-->"(" reached, end loop
          IF(br2.GT.0) GOTO 11

          !Procees only if char is not bracket -> valid line
          IF (br1.EQ.0) THEN

            j=j+1
            READ(OneLine,*) m%f(j)%owner

            !add +1 for difference in count start (Fortran starts with 1, OpenFOAM with 0)
            m%f(j)%owner = m%f(j)%owner + 1

          END IF

        END IF

        !Read first number in file -> defines num entries
        IF (isNumeric(KeyWord).AND.nentries.EQ.0) THEN
          READ(KeyWord,*) nentries
        END IF

      END DO
      
11    continue  
      IF (j.EQ.nentries) THEN
        PRINT *, "Foam owners read OK; nfaces=",j
      ELSE
        PRINT *, "Foam owners read nentries error; nfaces=",j
      END IF

      CLOSE(lun)

      CALL CPU_TIME(t2)
      PRINT *, "ReadFoamMesh elapsed time :: Owners:",t2-t1

! -----------------------------------------------------------------------------------------
!     Read face neighbours
! -----------------------------------------------------------------------------------------      

      lun = 13
      foamfile = "/neighbour"
      OPEN (lun,FILE=TRIM(fname) // TRIM(foamfile),ERR=10,STATUS='OLD')

      i=0
      j=0
      nentries=0

      OneLine="NOT EOF"
      DO WHILE (OneLine(1:3).NE.'EOF')

        i=i+1
        CALL rOneTL(lun,OneLine)
        READ (OneLine,*) KeyWord

        !PRINT *, "ol=",TRIM(OneLine)

        !If num entries > 0 (successfull number read) 
        IF(nentries.GT.0) THEN

          !Check if leading or terminating bracket is found
          br1=SCAN(TRIM(OneLine),"(")
          br2=SCAN(TRIM(OneLine),")")

          !When br2-->"(" reached, end loop
          IF(br2.GT.0) GOTO 12

          !Procees only if char is not bracket -> valid line
          IF (br1.EQ.0) THEN

            j=j+1
            READ(OneLine,*) m%f(j)%neighbour

            !add +1 for difference in count start (Fortran starts with 1, OpenFOAM with 0)
            IF (m%f(j)%neighbour.GE.0) THEN
              m%f(j)%neighbour = m%f(j)%neighbour + 1
            ELSE
              m%f(j)%neighbour = m%f(j)%neighbour
            END IF

          END IF

        END IF

        !Read first number in file -> defines num entries
        IF (isNumeric(KeyWord).AND.nentries.EQ.0) THEN
          READ(KeyWord,*) nentries
        END IF

      END DO

      12    continue  
!
!     When reading neighbour file, this nentries might be lower than num faces, important
!     is that j is equal to nentries value from neighbour file. If nentries < num faces, the bottom
!     correction is applied, see below.
!
      IF (j.EQ.nentries) THEN
        PRINT *, "Foam neighbours read OK; nfaces=",j
      ELSE
        PRINT *, "Foam neighbours read nentries error; nfaces=",j
      END IF
!
!     Foam mesh created with cfMesh marks boundary faces with -1 in neighbour file, to 
!     indicate no neighbours - > which means boundary face. With other meshes, in example
!     blockmesh, the boundary faces have no entries in neighbour file, and num entries
!     is smaller than num faces. If this is so (num entries < num faces), execute this loop
!     to fill -1 for boundary faces neighbour variable
! 
      IF ( j.LT.SIZE(m%f) ) THEN
        
        DO WHILE ( j.LT.SIZE(m%f) )
          
          j=j+1
          m%f(j)%neighbour = -1

        END DO

        PRINT *, "Foam neighbours adding -1 to neibs, now nfaces=",j

      END IF 

      CLOSE(lun)

      CALL CPU_TIME(t2)
      PRINT *, "ReadFoamMesh elapsed time :: Neighbours:",t2-t1
      !CALL BreakPoint()

      RETURN

10    continue ! error when opening input file
      CALL logWrite ("ERROR :: ReadFoamMesh :: Could not open : "//TRIM(fname)//TRIM(foamfile))
      CALL StopProgram(1)

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE ConstructFoamElementsFromFaces(m)
!     
!     Foam mesh is face based, use this routine to construct the mesh elements from 
!     faces
!        
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile   
         
      IMPLICIT NONE 

      TYPE(meshType) :: m

      INTEGER i,j,ilast,siz,ierr,maxElementFaces

      !ef -> elements from faces variable
      INTEGER, ALLOCATABLE :: ef(:,:)

      INTEGER k,l,iface,inode,nf,nn
      INTEGER, ALLOCATABLE :: nodecon(:)

      REAL t1,t2

      CALL CPU_TIME(t1)

      ierr = 0
      maxElementFaces = 100
!
!     ef(:,:) is an array of mesh elements, topologically 
!     connected with faces (kind of element-face connectivity)
!     Array structure:
!     ef(:,1) - last position in row, so we can 'append'
!     ef(:,2:maxElementFaces) - face indices
!      
      ALLOCATE(ef(m%nelem,maxElementFaces)) 
      
      ef(:,2:maxElementFaces) = 0

      !initialize first index with 2
      ef(:,1) = 2

      DO i = 1,m%nfac

        IF (m%f(i)%owner.GT.m%nelem) THEN
          ierr = 1
          GOTO 10
        ELSE IF (m%f(i)%neighbour.GT.m%nelem) THEN
          ierr = 2
          GOTO 10
        END IF

        !fill entries from owners (with i to indicate owner)
        ilast = ef(m%f(i)%owner,1)

        IF (ilast.GT.100) THEN
          PRINT *, "ConstructFoamElementsFromFaces :: More than 100 faces in element"
          PRINT *, "Stopping program, increase max elem size"
          CALL stopProgram(1)
        END IF

        ef( m%f(i)%owner,ilast ) = i
        ef( m%f(i)%owner,1 ) = ilast+1

        !fill entries from neighbours (with -i to indicate neighbour)
        IF (m%f(i)%neighbour.GT.-1) THEN
          ilast = ef(m%f(i)%neighbour,1)

          IF (ilast.GT.100) THEN
            PRINT *, "ConstructFoamElementsFromFaces :: More than 100 faces in element"
            PRINT *, "Stopping program, increase max elem size"
            CALL stopProgram(1)
          END IF

          ef( m%f(i)%neighbour,ilast ) = -i
          ef(m%f(i)%neighbour,1) = ilast+1   
        END IF

      END DO


      DO i = 1,m%nelem

        !ef(i,1) = size + 2, ilast from previous loop means last 
        !empty position, and index 1 is ilast value, so start is with 2
        siz = ef(i,1) - 2
        ALLOCATE(m%e(i)%faceCon(siz))
        m%e(i)%faceCon(:) = ef(i,2:siz+1)

      END DO

      !DO i = 1,20
      !  WRITE(*,'(A3,I4,A7,100I6)') "ef(",i,",1:10)=",&
      !  &( m%e(i)%faceCon(j), j=1,size(m%e(i)%faceCon) )
      !END DO

      DEALLOCATE(ef)

      !get element list of nodes (filter out repetitions)
      DO i = 1,m%nelem

        !number of faces in element
        nf = size( m%e(i)%faceCon )

        !init - max 100 nodes in element
        nn = 100

        ALLOCATE( nodecon(nn) )
        nodecon = -1

        l=0
        DO j = 1,nf

          iface = ABS( m%e(i)%faceCon(j) )
          DO k = 1,size( m%f(iface)%con )

            inode = m%f(iface)%con(k)
            IF ( .NOT. nodePresent(inode,nodecon,l) ) THEN
              l = l + 1
              nodecon(l) = inode
              !PRINT *, "Element ",i," has node ",inode
            END IF

          END DO

        END DO

        m%e(i)%nVertex = l

        ALLOCATE( m%e(i)%nodecon(m%e(i)%nVertex) )
        m%e(i)%nodecon = nodecon( 1:m%e(i)%nVertex )
        
        DEALLOCATE( nodecon )


      END DO


      CALL CPU_TIME(t2)
      PRINT *, "ConstructFoamElementsFromFaces elapsed time:",t2-t1

      RETURN

10    continue !Error
      PRINT *, "ERROR :: ConstructFoamElementsFromFaces :: ierr=",ierr
      CALL StopProgram(1)

      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      !check if node is already present in nodelist
      LOGICAL FUNCTION nodePresent(inode,nodelist,nn)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: inode, nn, nodelist(nn)
        INTEGER i

        nodePresent = .FALSE.

        IF ( nn.LT.1 ) RETURN

        DO i = 1,nn 

          IF ( inode.EQ.nodelist(i) ) THEN
            nodePresent = .TRUE.
            RETURN
          END IF

        END DO
    
      END FUNCTION nodePresent

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE TriangulateFoamFacesBackup(m)
!     
!     Divide faces read from foam file "faces" into triangles 
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile   
         
      IMPLICIT NONE 

      TYPE(meshType) :: m

      INTEGER i,j,l,ibest
      INTEGER np,ntri_local

      REAL(8) x1(3),x2(3),x3(3)

      INTEGER, ALLOCATABLE :: conn(:),tricon(:,:)
      REAL(8), ALLOCATABLE :: pnts(:,:)

      REAL t1,t2

      CALL CPU_TIME(t1)

      DO i = 1,m%nfac

        np = m%f(i)%nVertex
        ntri_local = np - 2

        ALLOCATE( m%f(i)%tf( ntri_local ) )

        ALLOCATE(conn(np),pnts(np,3),tricon(ntri_local,3))

        DO j=1,np
          conn(j) = m%f(i)%con(j)
          pnts(j,:) = m%x(conn(j),:)
        END DO

        IF (np .GT. 3) THEN

          !Calculate best angle to start decomposing
          !CalcAngle (array of point coordinates, number face vertex, return ibest)
          CALL calcAngles(pnts,np,ibest)
          CALL makeTriangles(pnts,conn,np,ibest,tricon)

        ELSE 

          tricon(ntri_local,:) = conn

        END IF

        !For all tri faces in a face
        DO j=1,ntri_local

          m%f( i )%tf( j )%con = tricon( j,: )

          ! ----------------------------------------------
          ! Get tri face center
          ! ----------------------------------------------

          DO l = 1,3

            x1(l) = m%x( m%f( i )%tf( j )%con(1),l )
            x2(l) = m%x( m%f( i )%tf( j )%con(2),l )
            x3(l) = m%x( m%f( i )%tf( j )%con(3),l )

            m%f( i )%tf( j )%center( l ) = ( 1.0D0/3.0D0 )*( x1(l) + x2(l) + x3(l)  )
    
          END DO

          ! ----------------------------------------------
          ! Get tri face normal
          ! ----------------------------------------------

          CALL CalNormal3p( m%x( m%f( i )%tf( j )%con(1),: ), &
          & m%x( m%f( i )%tf( j )%con(2),: ), &
          & m%x( m%f( i )%tf( j )%con(3),: ), &
          & m%f( i )%tf( j )%normal )

        END DO

        DEALLOCATE(conn,pnts,tricon)

      END DO

!      DO i = 1,m%nfac
!
!        PRINT *, "Face=",i
!        PRINT *, "nodes=",m%f(i)%nVertex
!        PRINT *, "numtf,",size( m%f( i )%tf )
!
!        DO j=1, size( m%f( i )%tf )
!          PRINT *, "con(",j,")=",m%f( i )%tf( j )%con
!        END DO
!
!        PRINT *, "-------------------"
!
!      END DO

      CALL CPU_TIME(t2)
      PRINT *, "TriangulateFoamFaces elapsed time:",t2-t1

      END   

! -----------------------------------------------------------------------------------------
      SUBROUTINE TriangulateFoamFaces(m)
!     
!     Divide faces read from foam file "faces" into triangles 
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile   
         
      IMPLICIT NONE 

      TYPE(meshType) :: m

      INTEGER i
      REAL t1,t2

      CALL CPU_TIME(t1)

      DO i = 1,m%nfac

        CALL doFaceTriangulation(m,i,0)

      END DO

!      DO i = 1,m%nfac
!
!        PRINT *, "Face=",i
!        PRINT *, "nodes=",m%f(i)%nVertex
!        PRINT *, "numtf,",size( m%f( i )%tf )
!
!        DO j=1, size( m%f( i )%tf )
!          PRINT *, "con(",j,")=",m%f( i )%tf( j )%con
!        END DO
!
!        PRINT *, "-------------------"
!
!      END DO

      CALL CPU_TIME(t2)
      PRINT *, "TriangulateFoamFaces elapsed time:",t2-t1

      END   

! -----------------------------------------------------------------------------------------
      SUBROUTINE doFaceTriangulation(m,iface,ibestoffset)
!     
!     Triangulate a iface
!     ibestoffset :: ibest = ibest + ibestoffset -> this should be more than 0
!     only when called from the correction routine, if it detects zerovolume geometry
!        
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile   
         
      IMPLICIT NONE 

      TYPE(meshType) :: m

      INTEGER i,j,l,ibest,iface,ibestoffset
      INTEGER np,ntri_local

      REAL(8) x1(3),x2(3),x3(3)

      INTEGER, ALLOCATABLE :: conn(:),tricon(:,:)
      REAL(8), ALLOCATABLE :: pnts(:,:)

      np = m%f(iface)%nVertex
      ntri_local = np - 2

      IF( .NOT. ASSOCIATED( m%f(iface)%tf )) THEN
        ALLOCATE( m%f(iface)%tf( ntri_local ) )
      END IF

      ALLOCATE(conn(np),pnts(np,3),tricon(ntri_local,3))

      DO j=1,np
        conn(j) = m%f(iface)%con(j)
        pnts(j,:) = m%x(conn(j),:)
      END DO

      IF (np .GT. 3) THEN

        !Calculate best angle to start decomposing
        !CalcAngle (array of point coordinates, number face vertex, return ibest)
        IF (ibestoffset.EQ.0) THEN
          CALL calcAngles(pnts,np,ibest)
        ELSE
          ibest = ibestoffset
        END IF

        CALL makeTriangles(conn,np,ibest,tricon)

      ELSE 

        tricon(ntri_local,:) = conn

      END IF

      !For all tri faces in a face
      DO j=1,ntri_local

        m%f( iface )%tf( j )%con = tricon( j,: )

        ! ----------------------------------------------
        ! Get tri face center
        ! ----------------------------------------------

        DO l = 1,3

          x1(l) = m%x( m%f( iface )%tf( j )%con(1),l )
          x2(l) = m%x( m%f( iface )%tf( j )%con(2),l )
          x3(l) = m%x( m%f( iface )%tf( j )%con(3),l )

          m%f( iface )%tf( j )%center( l ) = ( 1.0D0/3.0D0 )*( x1(l) + x2(l) + x3(l)  )
  
        END DO

        ! ----------------------------------------------
        ! Get tri face normal
        ! ----------------------------------------------

        CALL CalNormal3p( m%x( m%f( iface )%tf( j )%con(1),: ), &
        & m%x( m%f( iface )%tf( j )%con(2),: ), &
        & m%x( m%f( iface )%tf( j )%con(3),: ), &
        & m%f( iface )%tf( j )%normal )

      END DO

      DEALLOCATE(conn,pnts,tricon)


      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE CheckFoamTriangulation(m)
!
!     Find possible topological errors induced by triangulation
!
! -----------------------------------------------------------------------------------------

      USE mesh
          
      IMPLICIT NONE 

      TYPE(meshType) :: m

      INTEGER i
      INTEGER ierr
      INTEGER cnt, whcnt, unrepcnt

      cnt = 0
      unrepcnt = 0

      DO i = 1,m%nelem

        CALL CheckElementForZeroVolume(m,i,ierr)

        IF (ierr.NE.0) THEN

          cnt = cnt + 1

          whcnt = 0
          DO, WHILE (ierr.NE.0) 

            whcnt = whcnt + 1

            IF (ierr.NE.0) CALL RepairElementForZeroVolume(m,i,ierr) 
            CALL CheckElementForZeroVolume(m,i,ierr)
            
            IF (whcnt.GT.100) THEN
              unrepcnt = unrepcnt + 1
              PRINT *, "Element ",i," cannot repair zero volume"
              EXIT
            END IF

          END DO

        END IF

      END DO

      PRINT *, "Found ",cnt," elements with zerovol triangulation errors"
      PRINT *, "Successfully repaired ",cnt-unrepcnt
      PRINT *, "Remaining elements with zerovol triang error ",unrepcnt
      PRINT *, "Ending triangulation check"
      PRINT *, "----------------------------"

      IF (unrepcnt.GT.0) CALL BreakPoint()


      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE RepairElementForZeroVolume(m,ie,ierr)
!
! -----------------------------------------------------------------------------------------

      USE mesh
          
      IMPLICIT NONE 

      TYPE(meshType) :: m

      INTEGER ie,ierr
      INTEGER problemie,ie2
      INTEGER i,j,j2,k
      INTEGER f1,f2,iface,iface2
      INTEGER npairs,np,np2
      INTEGER ierrrep1,ierrrep2,ierrnei

      INTEGER stats(100,6)

      INTEGER, ALLOCATABLE :: cfaces(:,:) 

      npairs = size( m%conflictFaces(:,1) )
      ALLOCATE( cfaces(npairs,2) )

      cfaces = m%conflictFaces

      DO i = 1,npairs

        f1 = cfaces(i,1)
        f2 = cfaces(i,2)

        k = 0
        stats = 1

        iface = abs( m%e(ie)%faceCon(f1) )
        np = size( m%f(iface)%con )

        DO j = 1,np-1

          ierrrep1 = repairFace(ie,iface,j)

          iface2 = abs( m%e(ie)%faceCon(f2) )
          np2 = size( m%f(iface2)%con )

          DO j2 = 1,np2 - 1

            k = k + 1

            ierrrep2 = repairFace(ie,iface2,j2)

            !check nei elements status
            ierrnei = checkNeiElements(ie)

            stats(k,1) = ierrrep1
            stats(k,2) = ierrrep2
            stats(k,3) = ierrrep1 + ierrrep2
            stats(k,4) = size( m%inducedNeis )

            !successful repair
            IF (ierrrep1.EQ.0 .AND. ierrrep2.EQ.0 .AND. ierrnei.EQ.0) THEN
              GOTO 10
            END IF

          END DO
        END DO

        !DO j = 1,k
        !  PRINT *, "stats=",stats(j,:)
        !END DO

        !CALL WriteVTKTriangulatedElement("tricheck/trielem",1,m,ie)
        !CALL BreakPoint()

        !unsuccessful repair, return error
        ierr = 1
        !PRINT *, "Cannot repair, element=",ie," facepair=",i
        RETURN

10      CONTINUE !successfull repair, go check another pair if exists

      END DO

      PRINT *, "Successfully repaired element=",ie
      ierr = 0


      CONTAINS 

      INTEGER FUNCTION repairFace(ie,iface,j) 

        INTEGER,INTENT(IN) :: ie
        INTEGER,INTENT(IN) :: iface
        INTEGER,INTENT(IN) :: j

        INTEGER ierrown

        repairFace = 1

        !change triangulation order
        CALL doFaceTriangulation(m,iface,j)
    
        !check owner element status
        CALL CheckElementForZeroVolume(m,ie,ierrown)

        !successful repair
        IF (ierrown .EQ. 0) repairFace = 0

      END FUNCTION repairFace

      INTEGER FUNCTION checkNeiElements(ie)

        IMPLICIT NONE

        INTEGER,INTENT(IN) :: ie

        INTEGER i,iface,ienei,ierrnei
        INTEGER nnei

        INTEGER,ALLOCATABLE :: problemie(:)

        checkNeiElements = 0
        nnei = size(m%e(ie)%faceCon)

        ALLOCATE ( problemie(nnei) )
        problemie = 0

        DO i = 1,nnei

          iface = m%e(ie)%faceCon(i)

          IF (iface.GT.0) THEN
            ienei = m%f(abs(iface))%neighbour
          ELSE
            ienei = m%f(abs(iface))%owner
          END IF

          ierrnei = 0
          IF (ienei.GT.0) CALL CheckElementForZeroVolume(m,ienei,ierrnei)

          IF (ierrnei.GT.0) THEN 
            checkNeiElements = checkNeiElements + 1
            problemie(checkNeiElements) = ienei
          END IF

        END DO

        IF (ASSOCIATED(m%inducedNeis)) DEALLOCATE(m%inducedNeis)
        ALLOCATE( m%inducedNeis(checkNeiElements) )

        m%inducedNeis = problemie(1:checkNeiElements)


      END FUNCTION checkNeiElements

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE CheckElementForZeroVolume(m,ie,ierr)
!
!     Check triangles, so that 2 triangles of the same element do not
!     have the same 3 points. If they do, this is a zero volume shape problem
! -----------------------------------------------------------------------------------------

      USE mesh
          
      IMPLICIT NONE 

      TYPE(meshType) :: m

      INTEGER ie,ierr
      INTEGER i,iface 
      INTEGER nf,np,check

      INTEGER cf1,cf2 !conflict face pair local id's
      INTEGER, ALLOCATABLE :: cfaces(:,:)

      nf = size( m%e(ie)%faceCon )

      ALLOCATE( cfaces(nf,2) )

      ierr = 0
      DO i = 1,nf

        iface = abs( m%e(ie)%faceCon(i) )
        np = m%f(iface)%nVertex

        check = checkFace(ie,i,iface,cf1,cf2)

        IF (check.GT.0) THEN 
          ierr = ierr + 1
          cfaces(ierr,1) = cf1 
          cfaces(ierr,2) = cf2
        END IF

      END DO

      IF( ASSOCIATED(m%conflictFaces) ) DEALLOCATE(m%conflictFaces)
      ALLOCATE (m%conflictFaces(ierr,2))

      m%conflictFaces = cfaces(1:ierr,:)

      DEALLOCATE (cfaces)
      

      CONTAINS 
      INTEGER FUNCTION checkFace(ie,fj,iface,f1,f2)

        IMPLICIT NONE

        INTEGER,INTENT(IN) :: ie
        INTEGER,INTENT(IN) :: fj !local face id
        INTEGER,INTENT(IN) :: iface
        INTEGER,INTENT(OUT) :: f1 !conflict face local id 1
        INTEGER,INTENT(OUT) :: f2 !conflict face local id 2

        INTEGER :: k,j2,k2,iface2
        INTEGER :: n1,n2,n3,find1,find2,find3
        INTEGER :: FindLoc

        f1 = 0
        f2 = 0

        checkFace = 0

        DO k = 1,size( m%f(iface)%tf )

          !tri
          n1 = m%f(iface)%tf(k)%con(1)
          n2 = m%f(iface)%tf(k)%con(2)
          n3 = m%f(iface)%tf(k)%con(3)

          DO j2 = fj+1,size( m%e(ie)%faceCon )

            iface2 = abs( m%e(ie)%faceCon(j2) )
            DO k2 = 1,size( m%f(iface2)%tf )

              find1 = FindLoc( m%f(iface2)%tf(k2)%con, 3, n1 )
              find2 = FindLoc( m%f(iface2)%tf(k2)%con, 3, n2 )
              find3 = FindLoc( m%f(iface2)%tf(k2)%con, 3, n3 )

              IF( find1.NE.0 .AND. find2.NE.0 .AND. find3.NE.0 ) THEN

                !PRINT *, "Found problematic triangulation"
                !PRINT *, "Element ie=", ie
                !PRINT *, "Face iface1=",iface
                !PRINT *, "Tet itet1=",k
                !PRINT *, "Face iface2=",iface2
                !PRINT *, "Tet itet2=",k2
                !PRINT *, "n1,2,3=",n1,n2,n3
                !PRINT *, "con2=", m%f(iface2)%tf(k2)%con

                f1 = fj
                f2 = j2

                checkFace = 1
                RETURN           

              END IF

            END DO

          END DO

        END DO

      END FUNCTION checkFace

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE CalculateCentersNormals(m,mode)
!     
!     Calculate face / element center points and face normals
!        
!     Vol center mode:
!     mode = 1 :: Get vol center by averaging face center       
!     mode = 2 :: Get vol center by averaging element nodes
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile   
         
      IMPLICIT NONE 

      TYPE(meshType) :: m
      INTEGER mode

      INTEGER i,j,iface
      INTEGER npnts

      REAL(8),ALLOCATABLE :: x(:,:)
      REAL(8) xc(3),n(3)

      REAL t1,t2

      CALL CPU_TIME(t1)
!
!     Foam faces centers and normals
!
      DO i = 1,m%nfac
        
        npnts = size( m%f(i)%con )
        ALLOCATE(x( npnts,3 ))

        !get all face points into local pnts array
        DO j = 1,npnts
          x(j,:) = m%x( m%f(i)%con(j),: )
        END DO

        CALL CalPolygonCentroid(x,npnts,xc)
        m%f(i)%center = xc

        CALL CalPolygonNormal(x,npnts,n)
        m%f(i)%normal = n

        DEALLOCATE(x)

      END DO
!
!     Get foam elements volume centers
!
      IF (mode.EQ.1) THEN

        DO i = 1,m%nelem

          xc = 0.0D0
          DO j = 1,size( m%e(i)%faceCon )

            iface = ABS(m%e(i)%faceCon(j))
            xc = xc + m%f(iface)%center
            
          END DO

          m%e(i)%xc = xc / size( m%e(i)%faceCon )

        END DO

      END IF
!
!     Decomposed tets centers and normals
!
      DO i = 1,m%nelem

        DO j = 1,size( m%e(i)%tet )

          !PRINT *, "---------------------"
          !WRITE(*,'(A3,I6,A3,I6)') " i=",i," j=",j  
          !WRITE(*,'(A5,4I6)') "con=",m%e(i)%tet(j)%con
          !WRITE(*,'(A5,3F16.6)') "x1=",m%x(m%e(i)%tet(j)%con(1),:)
          !WRITE(*,'(A5,3F16.6)') "x2=",m%x(m%e(i)%tet(j)%con(2),:)
          !WRITE(*,'(A5,3F16.6)') "x3=",m%x(m%e(i)%tet(j)%con(3),:)
          !WRITE(*,'(A5,3F16.6)') "x4=",m%e(i)%xc

          ALLOCATE( x(3,3) )

          !WRITE(*,'(A5,4I6)') "con=",m%e(i)%tet(j)%con

          !------face 1-------!
          x(1,:) = m%x(  m%e(i)%tet(j)%con(2),: )
          x(2,:) = m%x(  m%e(i)%tet(j)%con(3),: )
          x(3,:) = m%e(i)%xc !point 4

          !WRITE(*,'(A3,9G12.4,2I6)') "x1=",x,i,j

          CALL CalPolygonCentroid(x,3,xc)
          m%e(i)%tet(j)%centers(1,:) = xc

          CALL CalPolygonNormal(x,3,n)
          m%e(i)%tet(j)%normals(1,:) = n

          !------face 2-------!
          x(1,:) = m%x(  m%e(i)%tet(j)%con(1),: )
          x(2,:) = m%e(i)%xc !point 4
          x(3,:) = m%x(  m%e(i)%tet(j)%con(3),: )

          !WRITE(*,'(A3,9G12.4,2I6)') "x2=",x,i,j

          CALL CalPolygonCentroid(x,3,xc)
          m%e(i)%tet(j)%centers(2,:) = xc

          CALL CalPolygonNormal(x,3,n)
          m%e(i)%tet(j)%normals(2,:) = n

          !------face 3-------!
          x(1,:) = m%x(  m%e(i)%tet(j)%con(1),: )
          x(2,:) = m%x(  m%e(i)%tet(j)%con(2),: )
          x(3,:) = m%e(i)%xc !point 4

          !WRITE(*,'(A3,9G12.4,2I6)') "x3=",x,i,j

          CALL CalPolygonCentroid(x,3,xc)
          m%e(i)%tet(j)%centers(3,:) = xc

          CALL CalPolygonNormal(x,3,n)
          m%e(i)%tet(j)%normals(3,:) = n

          !------face 4-------!
          x(1,:) = m%x(  m%e(i)%tet(j)%con(2),: )
          x(2,:) = m%x(  m%e(i)%tet(j)%con(1),: )
          x(3,:) = m%x(  m%e(i)%tet(j)%con(3),: )

          !WRITE(*,'(A3,9G12.4,2I6)') "x4=",x,i,j

          CALL CalPolygonCentroid(x,3,xc)
          m%e(i)%tet(j)%centers(4,:) = xc

          CALL CalPolygonNormal(x,3,n)
          m%e(i)%tet(j)%normals(4,:) = n

          DEALLOCATE(x)

        END DO

        !Check orientation of tet-normals
        CALL CheckTetNormals(m,i)

      END DO


      CALL CPU_TIME(t2)
      PRINT *, "CalculateCentersNormals elapsed time:",t2-t1

      END 


! -----------------------------------------------------------------------------------------
      SUBROUTINE CheckTetNormals(m,ie)
!     
!        
! -----------------------------------------------------------------------------------------
      
      USE mesh

      TYPE(meshType) :: m

      INTEGER ie,i,j,ntet
      REAL(8) c(3),n(3),p(3),res


      ntet = size ( m%e(ie)%tet )

      DO i = 1,ntet

        DO j = 1,4

          IF (j.LT.4) THEN
            p = m%x(  m%e(ie)%tet(i)%con(j),: )
          ELSE
            p = m%e(ie)%xc
          END IF

          c = m%e(ie)%tet(i)%centers(j,:)
          n = m%e(ie)%tet(i)%normals(j,:)

          !PRINT *, "p=",p
          !PRINT *, "c=",c
          !PRINT *, "n=",n

          CALL DotProduct(n,p-c,res)

          IF (res.GT.0.0D0) m%e(ie)%tet(i)%normals(j,:) = - n

        END DO

      END DO


      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE CreateFoamDecomposedTets(m)
!     
!     Create tets from face triangulation and element vol center         
!         
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile   
         
      IMPLICIT NONE 

      TYPE(meshType) :: m
      INTEGER i,j
      INTEGER iface,nf,ntet

      INTEGER,ALLOCATABLE :: tetcon(:,:)

      REAL t1,t2

      CALL CPU_TIME(t1)

      DO i = 1,m%nelem

        !count number of tetrahedrons - equal as element local nf*ntrif
        !
        ntet = 0
        nf = size( m%e(i)%faceCon )
        DO j = 1,nf
          iface = ABS( m%e(i)%faceCon(j) )
          ntet = ntet + size( m%f(iface)%tf )
        END DO

        !Build in-element tet connectivity
        !
        ALLOCATE ( tetcon(ntet,4) )
        tetcon = prepTets(i,ntet)

        ALLOCATE ( m%e(i)%tet(ntet) )
        DO j = 1,ntet
          m%e(i)%tet(j)%con = tetcon(j,:)
        END DO

        DEALLOCATE( tetcon )
        
      END DO

      CALL CPU_TIME(t2)
      PRINT *, "CreateFoamDecomposedTets elapsed time:",t2-t1


      ! --------------------------------
      CONTAINS    ! FUNCTIONS
      ! --------------------------------

      !
      ! Build in-element tet connectivity
      !

      FUNCTION prepTets(ie,ntet) RESULT (tetcon)

        IMPLICIT NONE

        INTEGER, INTENT(IN)   :: ie,ntet
        INTEGER i,j,k,nf,ntf,iface
        INTEGER tetcon(ntet,4)

        nf = size( m%e(ie)%faceCon )

        k=0
        DO i = 1,nf
          iface = ABS( m%e(ie)%faceCon(i) )

          ntf = size( m%f(iface)%tf )
          DO j = 1,ntf

            k=k+1

            tetcon(k,1:3) = m%f(iface)%tf(j)%con
            tetcon(k,4) = -1

          END DO

        END DO

      END FUNCTION



      END




! -----------------------------------------------------------------------------------------
      SUBROUTINE SetElementLC(m)
!
!     Calculate gradient of u
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure
!
!     Set local node coordinates
!

!
!           4               1 (0,0,0)
!          / \              2 (1,0,0) 
!         /   \   3         3 (0,1,0) 
!        /     \ /          4 (0,0,1)
!       1 ----- 2
!
      ALLOCATE(m%tetLC(4,3))
      m%tetLC(1,1)=0.0D0
      m%tetLC(1,2)=0.0D0
      m%tetLC(1,3)=0.0D0

      m%tetLC(2,1)=1.0D0
      m%tetLC(2,2)=0.0D0
      m%tetLC(2,3)=0.0D0

      m%tetLC(3,1)=0.0D0
      m%tetLC(3,2)=1.0D0
      m%tetLC(3,3)=0.0D0

      m%tetLC(4,1)=0.0D0
      m%tetLC(4,2)=0.0D0
      m%tetLC(4,3)=1.0D0      

!            zeta
!             ^
!             |                       -1 < xi,eta < 1
!             5                        0 < zeta < 1
!             |
!           1 ----- 4                  1 (-1,-1,0)
!          /  |    /                   2 (1, -1,0)
!         /   x - / ---> eta           3 (1,1,0)
!        /   /   /                     4 (-1,1,0)
!       2------ 3                      5 (0,0,1)
!          /
!         xi
      ALLOCATE(m%pyrLC(5,3))

      m%pyrLC(1,1)=-1.0D0
      m%pyrLC(1,2)=-1.0D0
      m%pyrLC(1,3)= 0.0D0

      m%pyrLC(2,1)= 1.0D0
      m%pyrLC(2,2)=-1.0D0
      m%pyrLC(2,3)= 0.0D0

      m%pyrLC(3,1)= 1.0D0
      m%pyrLC(3,2)= 1.0D0
      m%pyrLC(3,3)= 0.0D0

      m%pyrLC(4,1)=-1.0D0
      m%pyrLC(4,2)= 1.0D0
      m%pyrLC(4,3)= 0.0D0 

      m%pyrLC(5,1)= 0.0D0
      m%pyrLC(5,2)= 0.0D0
      m%pyrLC(5,3)= 1.0D0  

!           
!                 6         1 (0,0,-1)
!                /|         2 (1,0,-1)
!       4 ----- 5 |         3 (0,1,-1)
!       |       | |         4 (0,0,+1)
!       |       | 3         5 (1,0,+1)
!       |       |/          6 (0,1,+1)
!       1 ----- 2
!      
      ALLOCATE(m%priLC(6,3))

      m%priLC(1,1)= 0.0D0
      m%priLC(1,2)= 0.0D0
      m%priLC(1,3)=-1.0D0

      m%priLC(2,1)= 1.0D0
      m%priLC(2,2)= 0.0D0
      m%priLC(2,3)=-1.0D0

      m%priLC(3,1)= 0.0D0
      m%priLC(3,2)= 1.0D0
      m%priLC(3,3)=-1.0D0

      m%priLC(4,1)= 0.0D0
      m%priLC(4,2)= 0.0D0
      m%priLC(4,3)= 1.0D0 

      m%priLC(5,1)= 1.0D0
      m%priLC(5,2)= 0.0D0
      m%priLC(5,3)= 1.0D0  

      m%priLC(6,1)= 0.0D0
      m%priLC(6,2)= 1.0D0
      m%priLC(6,3)= 1.0D0  


!
!          8 ---- 7           1 (-1,-1,-1)
!         /|     /|
!        / |    / |
!       5 ---- 6  |
!       |  |   |  |
!       |  4 --|- 3
!       | /    |/
!       1 ---- 2
!      
      ALLOCATE(m%hexLC(8,3))

      m%hexLC(1,1)=-1.0D0
      m%hexLC(1,2)=-1.0D0
      m%hexLC(1,3)=-1.0D0

      m%hexLC(2,1)= 1.0D0
      m%hexLC(2,2)=-1.0D0
      m%hexLC(2,3)=-1.0D0

      m%hexLC(3,1)= 1.0D0
      m%hexLC(3,2)= 1.0D0
      m%hexLC(3,3)=-1.0D0

      m%hexLC(4,1)=-1.0D0
      m%hexLC(4,2)= 1.0D0
      m%hexLC(4,3)=-1.0D0 

      m%hexLC(5,1)=-1.0D0
      m%hexLC(5,2)=-1.0D0
      m%hexLC(5,3)= 1.0D0  

      m%hexLC(6,1)= 1.0D0
      m%hexLC(6,2)=-1.0D0
      m%hexLC(6,3)= 1.0D0  

      m%hexLC(7,1)= 1.0D0
      m%hexLC(7,2)= 1.0D0
      m%hexLC(7,3)= 1.0D0  

      m%hexLC(8,1)=-1.0D0
      m%hexLC(8,2)= 1.0D0
      m%hexLC(8,3)= 1.0D0  

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE SetGrad2(m,u,gradu)
!
!     Calculate gradient of u
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile  
      USE cpuTime    
      IMPLICIT NONE 


      TYPE(meshType) :: m  ! mesh data structure
      REAL(8) u(m%nnodes)

      REAL(8) gradu(m%nnodes,3) !final version, to be returned
      REAL(8) rawgradu(m%nnodes,4) !work version (:,1) is number of summations

      INTEGER i,j,inode
      REAL(8), ALLOCATABLE :: du(:,:)
      REAL(8) c ! cpu time
!
!     Init time measurement
!
      CALL cpStart(c)     
!
!     Init
!
      rawgradu = 0.0D0
      gradu=0.0D0
!
!     Loop over elements
!
      DO i=1,m%nelem

        ALLOCATE( du(m%e(i)%nVertex,3) )

        CALL SetElementGrad(m,u,i,du)  
        
        !add to gradOut matrix (summation for averaging later)
        DO j = 1,m%e(i)%nVertex

          inode = m%e(i)%nodecon(j)
          rawgradu(inode,1) = rawgradu(inode,1) + 1.0D0
          rawgradu(inode,2:4) = rawgradu(inode,2:4) + du(j,:)

        END DO  

        DEALLOCATE( du )

      END DO

      !average out all nodal gradients
      DO i = 1,m%nnodes
        gradu(i,1) = rawgradu(i,2) / rawgradu(i,1)
        gradu(i,2) = rawgradu(i,3) / rawgradu(i,1)
        gradu(i,3) = rawgradu(i,4) / rawgradu(i,1)
      END DO

      CALL cpStop(c)
      cput%meas(8)=cput%meas(8) + c


      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE SetElementGrad(m,u,ie,gradu)
!
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile  
      USE cpuTime    
      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure

      INTEGER ie,i,j

      REAL(8) u(m%nnodes)
      REAL(8) gradu(m%e(ie)%nVertex,3)

      INTEGER inode,iloc,ntet,FindLoc
      INTEGER tcon(4)

      REAL(8) gradOut(m%e(ie)%nVertex,4)
      REAL(8) uc,tx(4,3),tu(4),tetGradOut(4,3)

      !gradOut(:,1) - number of summations 
      !gradOut(:,2:4) - dx,dy,dz

      gradOut = 0.0D0
      gradu = 0.0D0

      !number of tetrahedrons
      ntet = size( m%e(ie)%tet )

      !get u in element center
      uc = 0.0D0
      DO i = 1,m%e(ie)%nVertex
        inode = m%e(ie)%nodeCon(i)
        uc = uc + u( inode )
      END DO
      uc = uc / m%e(ie)%nVertex

      !loop all tets
      DO i = 1,ntet

        tcon = m%e(ie)%tet(i)%con

        tx(1,:) = m%x( tcon(1),: )
        tx(2,:) = m%x( tcon(2),: )
        tx(3,:) = m%x( tcon(3),: )
        tx(4,:) = m%e(ie)%xc

        tu(1) = u( tcon(1) )
        tu(2) = u( tcon(2) )
        tu(3) = u( tcon(3) )
        tu(4) = uc

        CALL SimpleGrad(10,4,m,tx,tu,tetGradOut)

        !add to gradOut matrix (summation for averaging later)
        DO j = 1,3

          iloc = FindLoc( m%e(ie)%nodecon, m%e(ie)%nVertex, tcon(j) )

          gradOut(iloc,1) = gradOut(iloc,1) + 1.0D0
          gradOut(iloc,2:4) = gradOut(iloc,2:4) + tetGradOut(j,:)

        END DO  
        
      END DO

      !average gradOut matrix
      DO i = 1,m%e(ie)%nVertex
        gradu(i,1) = gradOut(i,2) / gradOut(i,1)
        gradu(i,2) = gradOut(i,3) / gradOut(i,1)
        gradu(i,3) = gradOut(i,4) / gradOut(i,1)
      END DO


      !PRINT *, "---------grad out-----------"
      !DO i = 1,m%e(ie)%nVertex
      !  WRITE( *,'(I3, 3F14.6)') INT(gradOut(i,1)), gradOut(i,2:4)
      !END DO

      END      
! -----------------------------------------------------------------------------------------
      SUBROUTINE SimpleGrad(type,nvert,m,x,u,gradOut)
!
!     Calculate gradient of u
!        
!     INPUTS:
!     type - vtk element type (tetrahedron is 10)
!     nvert - number of element vertices
!     m - mesh type
!     x(nvert,3) - element node coordinates matrix
!     u(nvert,3) - u values in element nodes
!
!     RETURNS:
!     gradOut(nvert,1) - du/dx
!     gradOut(nvert,2) - du/dy
!     gradOut(nvert,3) - du/dz
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile      
      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure

      INTEGER type,nvert   
      REAL(8) x(nvert,3),u(nvert)
      REAL(8) gradOut(nvert,3) 

      INTEGER i,ii,jj,kk
      REAL(8), ALLOCATABLE :: SFD(:,:),Jac(:,:),iJac(:,:)   
      REAL(8) dfildxyz

      ALLOCATE(Jac(3,3))  ! Jacobi matrix of derivatives
      ALLOCATE(iJac(3,3))  ! inverse Jacobi matrix of derivatives

!
!     Allocate SFD - shape function derivatives matrix
!
      ALLOCATE(SFD(nvert,3))
!
!     For all vertexes in an element
!
      DO i=1,nvert
!
!       global mesh node number
!      
        gradOut(i,:) = 0.0D0
!
!         Calculate shape function derivatives
!
        IF (type.EQ.12) THEN ! hexa  
          CALL hex_shapef_der(m%hexLC(i,:),SFD)                 
        ELSE IF (type.EQ.10) THEN ! tet
          CALL tet_shapef_der(m%tetLC(i,:),SFD)
        ELSE IF (type.EQ.14) THEN ! pyr
          CALL pyr_shapef_der(m%pyrLC(i,:),SFD)         
        ELSE IF (type.EQ.13) THEN ! prism
          CALL pri_shapef_der(m%priLC(i,:),SFD)          
        ELSE
          CALL logIntWrite ("ERROR :: SetGrad :: Element type not supported: ",type)
          CALL StopProgram(1)       
        END IF
!
!         Jacobi matrix of derivatives J_ij=\p (x,y,z) / \p (xi,eta,zeta)
!
        DO ii=1,3 ! Fx,Fy,Fz
          DO jj=1,3 ! xi,eta,zeta
            Jac(ii,jj)=0.0D0
            DO kk=1,nvert ! nodes
              !Jac(ii,jj) = Jac(ii,jj) + SFD(kk,jj) * m%x( conn(1,kk),ii )
              Jac(ii,jj) = Jac(ii,jj) + SFD(kk,jj) * x( kk,ii )
            END DO
          END DO 
        END DO
!
!         Calculate inverse of Jacobi matrix -> d(xi,eta,zeta)/d(x,y,z)
!
        CALL SetMatrixInverse3x3(Jac,iJac)
!
!         Calculate gradients for node 
!
        DO kk=1,nvert 
          DO ii=1,3 ! x,y,z              
            dfildxyz=0.0D0
            DO jj=1,3
              dfildxyz = dfildxyz + iJac(jj,ii)*SFD(kk,jj)
            END DO
            !gradOut(i,ii)=gradOut(i,ii)+u(conn(1,kk))*dfildxyz
            gradOut(i,ii)=gradOut(i,ii)+u( kk )*dfildxyz
          END DO
        END DO

      END DO ! nodes in element

      DEALLOCATE(SFD,Jac,iJac)

      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE ReadVTKfluid(m,fluid)
!
!     Read the flow information from the VTK file
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE mFluid
      USE logFile      
      IMPLICIT NONE 

      TYPE(meshType) :: m  ! mesh data structure
      TYPE(fluidType) :: fluid ! fluid results

      INTEGER lun

      lun = 12
      OPEN (lun,FILE=TRIM(m%FullName),ERR=10,STATUS='OLD')
      CALL logWrite ("MESSAGE :: ReadVTKfluid :: Reading : "//TRIM(m%FullName))

      CALL ReadVTKscalarField(lun,'CELL_DATA','p',m%nelem,fluid%Pe,fluid%iPe)
      CALL ReadVTKscalarField(lun,'POINT_DATA','p',m%nnodes,fluid%Pn,fluid%iPn)

      CALL ReadVTKvectorField(lun,'CELL_DATA','U',m%nelem,fluid%Ue,fluid%iUe)
      CALL ReadVTKvectorField(lun,'POINT_DATA','U',m%nnodes,fluid%Un,fluid%iUn)

      CALL ReadVTKvectorField(lun,'POINT_DATA','Vort',m%nnodes,fluid%Vortn,fluid%iVortn)
      CALL ReadVTKscalarField(lun,'POINT_DATA','T',m%nnodes,fluid%Tn,fluid%iTn)

      CLOSE (lun)
    
      RETURN

10    continue ! error when opening input file
      CALL logWrite ("ERROR :: ReadVTKfluid :: Could not open : "//TRIM(m%FullName))
      CALL StopProgram(1)
      
      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE ReadVTKvectorField(lun,s1,s2,n,f,isucc)
!
!     Reads a single scalar field from VTK file
!      
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER lun,n
      CHARACTER*(*) s1,s2
      REAL(8) f(n,3)
      CHARACTER KeyWord*64,OneLine*255,dummy*64      
      REAL(8), ALLOCATABLE :: tmp(:)
      INTEGER i,j,k,ipos,iSucc,iKey,ios
!     init
      f=0.0D0      
      iSucc=0
      iKey=0
!
!     Reading vector field
!
      REWIND(lun)
      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        READ (OneLine,*) KeyWord      
        IF (TRIM(KeyWord).EQ.TRIM(s1)) THEN
          iKey=1
          CALL rOneTL(lun,OneLine)
          DO WHILE (OneLine(1:3).NE.'EOF')
            READ (OneLine,*) KeyWord
            IF (TRIM(KeyWord).EQ.TRIM(s2)) THEN
              ALLOCATE (tmp(n*3))
              i=0
              DO WHILE (i<n*3)
                CALL rOneTL(lun,OneLine)  ! read whole line from file
11              i=i+1
                OneLine=ADJUSTL(OneLine)  ! remove leadin spaces
                !PRINT *, "ol/n=", i,n,3*n,OneLine
                READ(OneLine,*,IOSTAT=ios) tmp(i) ! read the first number in the line
                IF(ios.NE.0) THEN
                  IF(OneLine.EQ. 'EOF') THEN
                    GOTO 12
                  ELSE
                    PRINT *, "Significant error in ReadVTKvectorField Read"
                    CALL stopProgram(1)
                  ENDIF
                END IF
                ipos = SCAN(TRIM(OneLine)," ") ! find next space
                IF (ipos>0) THEN ! space found
                  OneLine=OneLine(ipos:LEN_TRIM(OneLine))
                  gOTO 11
                END IF
              END DO
12            k=0
              DO i=1,n
                DO j=1,3
                  k=k+1
                  f(i,j)=tmp(k) ! copy from temporary arry to mesh%x
                END DO
              END DO
              DEALLOCATE(tmp)
!             Field has been read
              iSucc=1              
            EXIT               
            END IF
            CALL rOneTL(lun,OneLine)
          END DO
        END IF
        IF (iKey.EQ.1) EXIT
        CALL rOneTL(lun,OneLine)
      END DO

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE ReadVTKscalarField(lun,s1,s2,n,f,iSucc)
!
!     Reads a single scalar field from VTK file
!      
! -----------------------------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER lun,n
      CHARACTER*(*) s1,s2
      REAL(8) f(n)
      CHARACTER KeyWord*64,OneLine*255,dummy*64      
      REAL(8), ALLOCATABLE :: tmp(:)
      INTEGER i,j,k,ipos,iSucc,iKey
!     init
      f=0.0D0      
      iSucc=0
      iKey=0
!
!     Reading CELL_DATA (scalar field)
!
      REWIND(lun)
      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        READ (OneLine,*) KeyWord
        IF (TRIM(KeyWord).EQ.TRIM(s1)) THEN
          iKey=1
          CALL rOneTL(lun,OneLine)
          DO WHILE (OneLine(1:3).NE.'EOF')
            READ (OneLine,*) KeyWord
            IF (TRIM(KeyWord).EQ.TRIM(s2)) THEN
              ALLOCATE (tmp(n))
              i=0
              DO WHILE (i<n)
                CALL rOneTL(lun,OneLine)  ! read whole line from file
11              i=i+1
                OneLine=ADJUSTL(OneLine)  ! remove leadin spaces
                READ(OneLine,*) tmp(i) ! read the first number in the line
                ipos = SCAN(TRIM(OneLine)," ") ! find next space
                IF (ipos>0) THEN ! space found
                  OneLine=OneLine(ipos:LEN_TRIM(OneLine))
                  gOTO 11
                END IF
              END DO
              k=0
              DO i=1,n
                k=k+1
                f(i)=tmp(k) ! copy from temporary arry to mesh%x
              END DO
              DEALLOCATE(tmp)
!             Field has been read
              iSucc=1
              EXIT               
            END IF
            CALL rOneTL(lun,OneLine)
          END DO
!         Point_Data keyword has been found
          EXIT
        END IF
        IF (iKey.EQ.1) EXIT        
        CALL rOneTL(lun,OneLine)
      END DO

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE ReadMeshPointsVTK(m,fname)
!
!     Read the mesh nodes from the VTK file      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile      
      IMPLICIT NONE 

      CHARACTER KeyWord*64,OneLine*255,dummy*64
      TYPE(meshType) :: m
      CHARACTER*(*) fname
      REAL(8), ALLOCATABLE :: tmp(:)
      INTEGER lun,i,j,k,ipos

      REAL t1,t2 

      CALL CPU_TIME(t1)

      m%FullName=fname

      lun = 12
      OPEN (lun,FILE=TRIM(m%FullName),ERR=10,STATUS='OLD')
!
!     Read
!
      CALL rOneTL(lun,OneLine)
      DO WHILE (OneLine(1:3).NE.'EOF')
        READ (OneLine,*) KeyWord
!
!       Reading POINTS
!
        IF (TRIM(KeyWord).EQ.'POINTS') THEN

          READ(OneLine,*) dummy,m%nnodes

          ALLOCATE (m%x(m%nnodes,3))
          !There are unknown number of datapoints in a sinle line
          ALLOCATE (tmp(m%nnodes*3))
          i=0
          DO WHILE (i<m%nnodes*3)
            CALL rOneTL(lun,OneLine)  ! read whole line from file
11          i=i+1
            OneLine=ADJUSTL(OneLine)  ! remove leadin spaces
            READ(OneLine,*) tmp(i) ! read the first number in the line
            ipos = SCAN(TRIM(OneLine)," ") ! find next space
            IF (ipos>0) THEN ! space found
              OneLine=OneLine(ipos:LEN_TRIM(OneLine))
              GOTO 11
            END IF
          END DO

          k=0
          DO i=1,m%nnodes
            DO j=1,3
              k=k+1
              m%x(i,j)=tmp(k) ! copy from temporary arry to mesh%x
            END DO
          END DO
          DEALLOCATE(tmp)
        END IF

!
!       Reading only number of elements, and allocating space
!
        IF (TRIM(KeyWord).EQ.'CELL_TYPES') THEN
          READ(OneLine,*) dummy,m%nelem
          ALLOCATE(m%e(m%nelem))
          GOTO 12
        END IF

        CALL rOneTL(lun,OneLine)
      END DO

12    CONTINUE

      CALL CPU_TIME(t2)
      PRINT *, "Read points VTK elapsed time:",t2-t1

      RETURN

10    CONTINUE ! error when opening input file
      CALL logWrite ("ERROR :: ReadVTKmesh :: Could not open : "//TRIM(m%FullName))
      CALL StopProgram(1)

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE WriteVTKflowResults(m,fluid)
!
!     Read the mesh information from the VTK file
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile      
      USE mFluid
      USE inpFile
      USE counters            

      IMPLICIT NONE 

      TYPE(meshType) :: m
      TYPE(fluidType) :: fluid ! fluid results
      CHARACTER(255) fname
      CHARACTER(255) vrstica,fmt
      INTEGER lun,i,j

!
!     Determine file name
!
      WRITE(fname,'(A,A,A,A,I0,A)') TRIM(inp%ParaviewResultsFolder),"/",TRIM(lg%IDname),".para-flow-",cnt%iTime,".vtk"
!
!     Output mesh part first
!
      CALL WriteVTKmesh(m,fname)
!
!     Append to file
!
      lun = 12
      OPEN (lun,FILE=TRIM(fname),ERR=10,status="old", position="append")

!
!     Export cell data
!      
      WRITE (lun,'(A,I10)') 'CELL_DATA ',m%nelem
      WRITE (lun,'(A,I10)') 'FIELD attributes 2'

      WRITE (lun,'(A,I10,A)') 'p 1 ',m%nelem," float"
      DO i=1,m%nelem
        WRITE (vrstica,*) REAL(fluid%Pe(i))
        CALL sqblnk(lun,vrstica)
      END DO      

      WRITE (lun,'(A,I10,A)') 'U 3 ',m%nelem," float"
      DO i=1,m%nelem
        WRITE (vrstica,*) REAL(fluid%Ue(i,1)),REAL(fluid%Ue(i,2)),REAL(fluid%Ue(i,3))
        CALL sqblnk(lun,vrstica)
      END DO      
!
!     Export field data
!      
      WRITE (lun,'(A,I10)') 'POINT_DATA ',m%nnodes
      WRITE (lun,'(A,I10)') 'FIELD attributes 2'

      WRITE (lun,'(A,I10,A)') 'p 1 ',m%nnodes," float"
      DO i=1,m%nnodes
        WRITE (vrstica,*) REAL(fluid%Pn(i))
        CALL sqblnk(lun,vrstica)
      END DO      

      WRITE (lun,'(A,I10,A)') 'U 3 ',m%nnodes," float"
      DO i=1,m%nnodes
        WRITE (vrstica,*) REAL(fluid%Un(i,1)),REAL(fluid%Un(i,2)),REAL(fluid%Un(i,3))
        CALL sqblnk(lun,vrstica)
      END DO      

!      WRITE (lun,'(A,I10,A)') 'Vort 3 ',m%nnodes," float"
!      DO i=1,m%nnodes
!        WRITE (vrstica,*) REAL(fluid%Vortn(i,1)),REAL(fluid%Vortn(i,2)),REAL(fluid%Vortn(i,3))
!        CALL sqblnk(lun,vrstica)
!      END DO      

!      WRITE (lun,'(A,I10,A)') 'T 1 ',m%nnodes," float"
!      DO i=1,m%nnodes
!        WRITE (vrstica,*) REAL(fluid%Tn(i))
!        CALL sqblnk(lun,vrstica)
!      END DO      


      CLOSE (lun)

      RETURN


10    continue ! error when opening input file
      CALL logWrite ("ERROR :: WriteVTKflowResults :: Could not open : "//TRIM(fname))
      CALL StopProgram(1)
      
      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE WriteVTKmesh(m,fname)
!
!     Read the mesh information from the VTK file
!      
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile      
      IMPLICIT NONE 

      TYPE(meshType) :: m
      CHARACTER*(*) fname
      CHARACTER(255) vrstica,fmt
      INTEGER lun,i,j

      lun = 12

      !This will work for only non-polyhedra
      !And even this will probably not work, since we are not using this approach anymore
      !Proably best to rewrie in terms of polyhedrons

      OPEN (lun,FILE=TRIM(fname),ERR=10,STATUS='UNKNOWN')
      WRITE (lun,'(A)') '# vtk DataFile Version 2.0'
      WRITE (lun,'(A)') 'SuperV VTK mesh export'
      WRITE (lun,'(A)') 'ASCII'
      WRITE (lun,'(A)') 'DATASET UNSTRUCTURED_GRID'
!
!     Export points
!      
      WRITE (lun,'(A,I10,A)') 'POINTS',m%nnodes," float"
      DO i=1,m%nnodes
        WRITE (vrstica,*) REAL(m%x(i,1)),REAL(m%x(i,2)),REAL(m%x(i,3))
        CALL sqblnk(lun,vrstica)
      END DO
!
!     Export cell conectivity
!      
      WRITE (lun,'(A,I10,I10)') 'CELLS ',m%nelem,m%ncon  
      DO i=1,m%nelem
        WRITE (fmt,'(A,I2,A)') "(I10,",m%e(i)%nVertex,"(I10))"
        WRITE (vrstica,fmt) m%e(i)%nVertex,(m%e(i)%con(1,j)-1,j=1,m%e(i)%nVertex) ! +1 : zero start con list
        CALL sqblnk(lun,vrstica)
      END DO
!
!     Export cell types
!      
      WRITE (lun,'(A,I10)') 'CELL_TYPES ',m%nelem
      DO i=1,m%nelem
        WRITE (vrstica,fmt) m%e(i)%type
        CALL sqblnk(lun,vrstica)
      END DO

      CLOSE (lun)

      RETURN


10    continue ! error when opening input file
      CALL logWrite ("ERROR :: WriteVTKmesh :: Could not open : "//TRIM(fname))
      CALL StopProgram(1)
      
      END  

! -----------------------------------------------------------------------------------------
      SUBROUTINE WriteVTKTriangle(fname,x1,x2,x3,center,normal)
!
!     Read the mesh information from the VTK file
!      
! -----------------------------------------------------------------------------------------
  
      USE logFile  
      IMPLICIT NONE       

      CHARACTER*(*) fname
      INTEGER lun
      REAL(8) x1(3),x2(3),x3(3),center(3),normal(3)

      REAL(8) d,d1,d2,d3

      lun = 12

      !get average edge, for normal vector scaling
      CALL dist2P(x2,x3,d1)
      CALL dist2P(x3,x1,d2)
      CALL dist2P(x1,x2,d3)

      d = (1.0D0/3.0D0)*(d1+d2+d3) 


      OPEN (lun,FILE=TRIM(fname),ERR=10,STATUS='UNKNOWN')
      WRITE (lun,'(A)') '# vtk DataFile Version 2.0'
      WRITE (lun,'(A)') 'SuperV VTK mesh export'
      WRITE (lun,'(A)') 'ASCII'
      WRITE (lun,'(A)') 'DATASET UNSTRUCTURED_GRID'

      WRITE (lun,'(A,I10,A)') 'POINTS ',5," float"
      WRITE (lun,'(3F14.6)') x1
      WRITE (lun,'(3F14.6)') x2
      WRITE (lun,'(3F14.6)') x3
      WRITE (lun,'(3F14.6)') center
      WRITE (lun,'(3F14.6)') center+0.25*d*normal

      WRITE (lun,'(A,I10,I10)') 'CELLS ', 2, 7
      WRITE (lun,'(A)') '3 0 1 2'
      WRITE (lun,'(A)') '2 3 4'

      WRITE (lun,'(A,I10)') 'CELL_TYPES ', 2
      WRITE (lun,'(A)') '5'
      WRITE (lun,'(A)') '3'

      CLOSE (lun)

      RETURN


10    continue ! error when opening input file
      CALL logWrite ("ERROR :: WriteVTKTriangle :: Could not open : "//TRIM(fname))
      CALL StopProgram(1)
      
      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE WriteVTKPoint(fname,iT,x)
!
!     Read the mesh information from the VTK file
!      
! -----------------------------------------------------------------------------------------
  
      USE logFile  
      IMPLICIT NONE       

      CHARACTER*(*) fname
      CHARACTER*1024 filename
      CHARACTER*20 itstring

      INTEGER lun,iT
      REAL(8) x(3)

      lun = 12

      WRITE(itstring,'(I10)') iT 
      WRITE(filename,'(A,A,A)') TRIM(fname),TRIM(ADJUSTL(itstring)),'.vtk'

      OPEN (lun,FILE=TRIM(filename),ERR=10,STATUS='UNKNOWN')
      WRITE (lun,'(A)') '# vtk DataFile Version 2.0'
      WRITE (lun,'(A)') 'SuperV VTK mesh export'
      WRITE (lun,'(A)') 'ASCII'
      WRITE (lun,'(A)') 'DATASET UNSTRUCTURED_GRID'

      WRITE (lun,'(A,I10,A)') 'POINTS ',1," float"
      WRITE (lun,'(3F14.6)') x

      WRITE (lun,'(A,I10,I10)') 'CELLS ', 1, 2
      WRITE (lun,'(A)') '1 0'

      WRITE (lun,'(A,I10)') 'CELL_TYPES ', 1
      WRITE (lun,'(A)') '1'


      CLOSE (lun)

      RETURN


10    continue ! error when opening input file
      CALL logWrite ("ERROR :: WriteVTKPoint :: Could not open : "//TRIM(fname))
      CALL StopProgram(1)
      
      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE WriteVTKTriangulatedFace(fname,iT,m,iface)
!
!      
! -----------------------------------------------------------------------------------------
  
      USE logFile 
      USE mesh 

      IMPLICIT NONE  
      TYPE(meshType) :: m  ! mesh data structure

      CHARACTER*(*) fname
      CHARACTER*1024 filename
      CHARACTER*20 itstring

      INTEGER ie, lun, np, iface, iT
      INTEGER FindLoc
      REAL(8) x1(3),x2(3),x3(3),center(3),normal(3)

      REAL(8) d,d1,d2,d3

      INTEGER i,j,k

      INTEGER, ALLOCATABLE :: lcon(:)
      INTEGER :: tricon(4)

      lun = 12

      !Write header

      WRITE(itstring,'(I10)') iT 
      WRITE(filename,'(A,A,A)') TRIM(fname),TRIM(ADJUSTL(itstring)),'.vtk'

      OPEN (lun,FILE=TRIM(filename),ERR=10,STATUS='UNKNOWN')
      WRITE (lun,'(A)') '# vtk DataFile Version 2.0'
      WRITE (lun,'(A)') 'SuperV VTK mesh export'
      WRITE (lun,'(A)') 'ASCII'
      WRITE (lun,'(A)') 'DATASET UNSTRUCTURED_GRID'

      IF (ie.LT.1) THEN
        WRITE (lun,'(A)') '---------------'
        WRITE (lun,'(A)') 'Error iface < 0'
        WRITE (lun,'(A)') '---------------'
        GOTO 11
      END IF

      !Write points
      np = size( m%f(iface)%con )

      ALLOCATE( lcon(np) )
      WRITE (lun,'(A,I6,A)') 'POINTS ',np," float"
      DO i = 1,np
        lcon(i) = m%f(iface)%con(i)
        WRITE (lun,'(3F14.6)') m%x( lcon(i),: )
      END DO

      !Write triangle cells
      k = size( m%f(iface)%tf )
      WRITE (lun,'(A,I6,I6)') 'CELLS ', k, 4*k

      DO i = 1,k 

          tricon(1) = 3
          tricon(2) = FindLoc( lcon, np, m%f(iface)%tf(i)%con(1) ) - 1
          tricon(3) = FindLoc( lcon, np, m%f(iface)%tf(i)%con(2) ) - 1
          tricon(4) = FindLoc( lcon, np, m%f(iface)%tf(i)%con(3) ) - 1
          WRITE (lun,'(4I6)') tricon

      END DO

      WRITE (lun,'(A,I10)') 'CELL_TYPES ', k

      DO i = 1,k
        WRITE (lun,'(A)') '5'
      END DO 

11    CONTINUE     

      CLOSE (lun)
      RETURN


10    continue ! error when opening input file
      CALL logWrite ("ERROR :: WriteVTKTriangulatedFace:: Could not open : "//TRIM(fname))
      CALL StopProgram(1)
      
      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE WriteVTKTriangulatedElement(fname,iT,m,ie)
!
!     Read the mesh information from the VTK file
!      
! -----------------------------------------------------------------------------------------
  
      USE logFile 
      USE mesh 

      IMPLICIT NONE  
      TYPE(meshType) :: m  ! mesh data structure

      CHARACTER*(*) fname
      CHARACTER*1024 filename
      CHARACTER*20 itstring

      INTEGER ie, lun, np, iface, iT
      INTEGER FindLoc
      REAL(8) x1(3),x2(3),x3(3),center(3),normal(3)

      REAL(8) d,d1,d2,d3

      INTEGER i,j,k

      INTEGER, ALLOCATABLE :: lcon(:)
      INTEGER :: tricon(4)

      lun = 12

      !Write header

      WRITE(itstring,'(I10)') iT 
      WRITE(filename,'(A,A,A)') TRIM(fname),TRIM(ADJUSTL(itstring)),'.vtk'

      OPEN (lun,FILE=TRIM(filename),ERR=10,STATUS='UNKNOWN')
      WRITE (lun,'(A)') '# vtk DataFile Version 2.0'
      WRITE (lun,'(A)') 'SuperV VTK mesh export'
      WRITE (lun,'(A)') 'ASCII'
      WRITE (lun,'(A)') 'DATASET UNSTRUCTURED_GRID'

      IF (ie.LT.1) THEN
        WRITE (lun,'(A)') '---------------'
        WRITE (lun,'(A)') 'Error ie < 0'
        WRITE (lun,'(A)') '---------------'
        GOTO 11
      END IF

      !Write points
      np = size( m%e(ie)%nodeCon )

      ALLOCATE( lcon(np) )
      WRITE (lun,'(A,I6,A)') 'POINTS ',np," float"
      DO i = 1,np
        lcon(i) = m%e(ie)%nodeCon(i)
        WRITE (lun,'(3F14.6)') m%x( lcon(i),: )
      END DO

      !Write triangle cells

      k = 0
      DO i = 1,size(m%e(ie)%faceCon)

        iface = abs(m%e(ie)%faceCon(i))
        k = k + size(m%f(iface)%tf)

      END DO

      WRITE (lun,'(A,I6,I6)') 'CELLS ', k, 4*k

      DO i = 1,size( m%e(ie)%faceCon )
        iface = abs(m%e(ie)%faceCon(i))
        DO j = 1,size( m%f(iface)%tf )

          tricon(1) = 3
          tricon(2) = FindLoc( lcon, np, m%f(iface)%tf(j)%con(1) ) - 1
          tricon(3) = FindLoc( lcon, np, m%f(iface)%tf(j)%con(2) ) - 1
          tricon(4) = FindLoc( lcon, np, m%f(iface)%tf(j)%con(3) ) - 1
          WRITE (lun,'(4I6)') tricon

        END DO
      END DO

      WRITE (lun,'(A,I10)') 'CELL_TYPES ', k

      DO i = 1,k
        WRITE (lun,'(A)') '5'
      END DO 

11    CONTINUE     

      CLOSE (lun)
      RETURN


10    continue ! error when opening input file
      CALL logWrite ("ERROR :: WriteVTKTriangulatedElement :: Could not open : "//TRIM(fname))
      CALL StopProgram(1)
      
      END