!
! -----------------------------------------------------------------------------------------
!
      MODULE mesh
!
! -----------------------------------------------------------------------------------------
!     
      !For now this is sphere mesh       
      TYPE localMeshType
        CHARACTER*255 ::  FullName            ! path & file name
        INTEGER :: nbelem,nicell,nnodes,nbnodes,npob
        
        INTEGER, POINTER :: ibc(:,:)  ! nbelem,npob
        REAL(8), POINTER :: x(:,:)    ! nnodes,npx

      END TYPE localMeshType 
!
! -----------------------------------------------------------------------------------------
!
      TYPE faceType

        INTEGER :: nVertex ! number of nodes in element
        INTEGER :: owner !face owner element
        INTEGER :: neighbour !face neighbour element
        INTEGER, POINTER :: con(:) ! face connectivity (nVertex)
        REAL(8) :: center(3) !face center
        REAL(8) :: normal(3) !face normals

        TYPE(triFaceType), POINTER :: tf(:) !type for triangles (face decomposition)

      END TYPE faceType
!
! -----------------------------------------------------------------------------------------
!
      TYPE triFaceType

        INTEGER :: con(3) ! face connectivity (nVertex)
        REAL(8) :: center(3) !face center
        REAL(8) :: normal(3) !face normals

      END TYPE triFaceType
!
! -----------------------------------------------------------------------------------------
!      

      TYPE elementType
        INTEGER :: type   ! type of element 
        INTEGER :: nVertex ! number of nodes in element
        !INTEGER :: nNeib  ! number of neihourin elements
        INTEGER, POINTER :: con(:,:) ! connectivity (nVertex)
        !INTEGER, POINTER :: neib(:) ! list of neibours (nNeib)  
        REAL(8) :: xc(3) !element vol center, from foam data
!
!       Connectivities
! 
        INTEGER, POINTER  :: faceCon(:) ! element faces list
        INTEGER, POINTER  :: nodeCon(:) ! list of unique vertex ID, for polyhedron

        TYPE(tetElementType), POINTER :: tet(:) !type for in-element tetrahedrons

      END TYPE elementType 
!
! -----------------------------------------------------------------------------------------
! 
      TYPE tetElementType

        INTEGER :: con(4) ! tet connectivity
        REAL(8) :: centers(4,3)
        REAL(8) :: normals(4,3)

      END TYPE tetElementType     
!
! -----------------------------------------------------------------------------------------
! 
      TYPE boundaryType

        INTEGER :: type       ! 0 - default 
                              ! 1 - wall
                              ! 2 - patch
        CHARACTER*20 :: name
        INTEGER :: startFace
        INTEGER :: nFaces

      END TYPE boundaryType 
!
! -----------------------------------------------------------------------------------------
!      
      TYPE meshType

        CHARACTER*255 ::  FullName            ! path & file name
        
        INTEGER :: nelem  ! number of mesh elements
        INTEGER :: nnodes ! number of nodes
        INTEGER :: ncon   ! number of inteers in cell connectivity
        INTEGER :: nbound ! number of mesh boundaries
        REAL(8) xmin,xmax,ymin,ymax,zmin,zmax ! mesh extents
        
        TYPE(elementType), POINTER :: e(:)  ! list of elements (nelem) 
        REAL(8), POINTER :: x(:,:)    ! nnodes,3

!       Local coordinates of nodes in mesh element
        REAL(8), POINTER :: tetLC(:,:)    ! 4,3 (4nodes,xez)
        REAL(8), POINTER :: priLC(:,:)    ! 6,3 (6nodes,xez)
        REAL(8), POINTER :: pyrLC(:,:)    ! 5,3 (5nodes,xez)
        REAL(8), POINTER :: hexLC(:,:)    ! 8,3 (8nodes,xez)
!
!       Fast locate particle cell 
!        
        INTEGER, POINTER  :: boxElemList(:) => null()!list of initial search elements (faster find particles at t=0)
        REAL(8), POINTER  :: lpc_w(:) !nelem ! vsota koordinat        
        INTEGER, POINTER  :: lpc_sortw(:) ! sortirani ID-ji elementov na velikost w
        REAL(8)           :: lpc_wrmax ! maksimalni w radij celotne mreze=3*rmax?    
!
!       Foam mesh variables
! 
        INTEGER :: nfac   ! number of mesh faces        
        TYPE(faceType), POINTER :: f(:) ! list of faces (nfac) 
        TYPE(triFaceType), POINTER :: tf(:) ! list of decomposed triangular faces (ntfac)
        TYPE(boundaryType), POINTER :: bound(:) !list of boundaries
!
!       Foam check triangulation variables
!         
        INTEGER,POINTER :: conflictFaces(:,:)
        INTEGER,POINTER :: inducedNeis(:)
        
      END TYPE meshType 

      END MODULE