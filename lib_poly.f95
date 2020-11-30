! -----------------------------------------------------------------------------------------
      SUBROUTINE GetPolyTetDecomWeights(m,elid,v,lambda_i,tcon,ierr)
!
!     Mean value coordinates method for in polyhedron interpolation
!     m --> TYPE(VTKmeshType)
!     elid --> element ID [m%e(elid)]
!     v --> point in which to interpolate nodal values
!     
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile      
      IMPLICIT NONE 
      TYPE(VTKmeshType) :: m
  
      INTEGER elid,tcon(3)
      REAL(8) v(3)

      INTEGER i,j,k1,k2,ifound,jmin
      INTEGER ntri,nv,found,ierr

      REAL(8) res,dist,mindist
      REAL(8) xez(3,3),xtri(4,3)
      !REAL(8) lambda_i(4,3)
      REAL(8) lambda_i(4)

      found = 0
      ntri=size(m%e(elid)%triCon(:,1))
      nv=m%e(elid)%nVertex

!-----------------------------------------------------------      
!
!       Find in which tet
!
!----------------------------------------------------------- 

      !PRINT *, "elid=",elid
      !WRITE(*,'(A2,3F10.5)') "v=",v

      DO i = 1,ntri
        !PRINT *, "---------------------------------------------------------"
        !WRITE(*,'(A4,4I6)') "con=",m%e(elid)%tricon(i,:)
        !WRITE(*,'(A4,3F10.5)') "x1=",m%x( m%e(elid)%tricon(i,1) ,:)
        !WRITE(*,'(A4,3F10.5)') "x2=",m%x( m%e(elid)%tricon(i,2) ,:)
        !WRITE(*,'(A4,3F10.5)') "x3=",m%x( m%e(elid)%tricon(i,3) ,:)
        !WRITE(*,'(A4,3F10.5)') "x4=",m%e(elid)%volCenter
        !PRINT *, ""

        DO j = 1,4
          
          k1 = 3*(j-1)+1
          k2 = 3*j
          
          CALL DotProduct( m%e(elid)%tetN(i,k1:k2), v - m%e(elid)%tetC(i,k1:k2), res )
          !WRITE(*,'(A4,I5,A4,I5,A4,I5)') "  i=",i,"  k1=",k1,"  k2=",k2
          !WRITE(*,'(A4,3F10.5)') "tetC=",m%e(elid)%tetC(i,k1:k2)
          !WRITE(*,'(A4,3F10.5)') "tetN=",m%e(elid)%tetN(i,k1:k2)
          !WRITE(*,'(A4,3F10.5)') "parR=",v - m%e(elid)%tetC(i,k1:k2)
          !WRITE(*,'(A4,F10.5)') "res=",res
          
          IF (res.GT.0.0D0) THEN
            GOTO 12
          END IF
        END DO
        found = found + 1
        ifound = i

12      res=res
      END DO

      IF (found.EQ.1) THEN
        !PRINT *, "TetDecom:: Message:: Particle found once, all good"
        m%iff = m%iff + 1
      ELSE IF (found.EQ.0) THEN

        ierr = 2

        !PRINT *, "TetDecom:: Error:: Particle not found at Elid:", elid
        m%inf = m%inf + 1

        !WRITE(*,'(A2,3F10.3)') "x=",m%x ( m%e(elid)%tricon(i,1) ,:)
        !WRITE(*,'(A2,3F10.3)') "x=",m%x ( m%e(elid)%tricon(i,2) ,:)
        !WRITE(*,'(A2,3F10.3)') "x=",m%x ( m%e(elid)%tricon(i,3) ,:)
        !WRITE(*,'(A2,3F10.3)') "x=",m%e(elid)%volCenter

        mindist = 9.99E+10
        DO j = 1,m%e(elid)%nVertex

          CALL vecLen(v-m%x ( m%e(elid)%polyUniver(j) ,:),dist)
          !PRINT *, "dist=",dist
          IF (dist.LT.mindist) THEN
            mindist = dist 
            jmin = j
          END IF        
          !WRITE(*,'(A2,3F12.6)') "x=",m%x ( m%e(elid)%polyUniver(j) ,:)

        END DO

        !PRINT *, "mindist=",mindist
        !WRITE(*,'(A2,3F12.6)') "pnt=",v
        !CALL stopProgram(1)

        GOTO 10

      ELSE IF (found.GT.1) THEN

        ierr = 3
        !PRINT *, "TetDecom:: Error:: Elid:",elid," Particle multiple found, n=",found
        m%imf = m%imf + 1
        !CALL StopProgram(1)
      ELSE
        PRINT *, "TetDecom:: Fatal Error:: Impossible error, go home fortran"
        CALL StopProgram(1)
      END IF

!-----------------------------------------------------------      
!       Get tet con, xez and shape functions
!----------------------------------------------------------- 

      tcon(1) = m%e(elid)%triCon(ifound,1)
      tcon(2) = m%e(elid)%triCon(ifound,2)
      tcon(3) = m%e(elid)%triCon(ifound,3)

      xtri(1,:) = m%x( m%e(elid)%triCon(ifound,1) ,:)
      xtri(2,:) = m%x( m%e(elid)%triCon(ifound,2) ,:)
      xtri(3,:) = m%x( m%e(elid)%triCon(ifound,3) ,:)
      xtri(4,:) = m%e(elid)%volCenter(:)
      CALL tet_kks2lks(v,xtri(1,:),xtri(2,:),xtri(3,:),xtri(4,:),xez,ierr)

      !WRITE(*,'(3F9.2)') xez(1,:)
      !WRITE(*,'(3F9.2)') xez(2,:)
      !WRITE(*,'(3F9.2)') xez(3,:)

      CALL tet_shapef(xez,lambda_i)

      RETURN

10    CONTINUE

      tcon(1) = jmin
      tcon(2) = m%e(elid)%triCon(1,2)
      tcon(3) = m%e(elid)%triCon(1,3)

      lambda_i = 0.0D0
      lambda_i(1) = 1.0D0

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE CalcPolyTetCentNor(m)
!
!     Calculate "in-poly" tet centers and normals
!     
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile      
      IMPLICIT NONE 
      TYPE(VTKmeshType) :: m
  
      INTEGER i,j,k
      INTEGER ntri,npoly

      REAL(8) xtri(4,3)

      npoly = 0
      DO i=1,m%nelem     
      !DO i = 1329,1329   
        IF (m%e(i)%type.EQ.42) THEN

          ntri = size(m%e(i)%triCon(:,1))
          ALLOCATE(m%e(i)%tetC(ntri,12))
          ALLOCATE(m%e(i)%tetN(ntri,12))
!
!         For all triangles
!          
          DO j=1,ntri
          !DO j=10,10

            !PRINT *, "--------------------------------------------------"
            !PRINT *, "j=",j

            xtri(1,:) = m%x( m%e(i)%triCon(j,1) ,:)
            xtri(2,:) = m%x( m%e(i)%triCon(j,2) ,:)
            xtri(3,:) = m%x( m%e(i)%triCon(j,3) ,:)
            xtri(4,:) = m%e(i)%volCenter(:)
            !WRITE(*,'(A4,3F10.5)') "x1=",xtri(1,:)
            !WRITE(*,'(A4,3F10.5)') "x2=",xtri(2,:)
            !WRITE(*,'(A4,3F10.5)') "x3=",xtri(3,:)
            !WRITE(*,'(A4,3F10.5)') "x4=",xtri(4,:)
!
!           Get tet centers
!
            DO k=1,3
              m%e(i)%tetC(j,6+k) = (1.0D0/3.0D0) * ( xtri(1,k) + xtri(2,k) + xtri(4,k) ) !3!
              m%e(i)%tetC(j,0+k) = (1.0D0/3.0D0) * ( xtri(2,k) + xtri(3,k) + xtri(4,k) ) !1!
              m%e(i)%tetC(j,3+k) = (1.0D0/3.0D0) * ( xtri(1,k) + xtri(4,k) + xtri(3,k) ) !2!
              m%e(i)%tetC(j,9+k) = (1.0D0/3.0D0) * ( xtri(2,k) + xtri(1,k) + xtri(3,k) ) !4!
            END DO

            !WRITE(*,'(A4,3F10.5)') "c1=",m%e(i)%tetC(j,1:3)
            !WRITE(*,'(A4,3F10.5)') "c2=",m%e(i)%tetC(j,4:6)
            !WRITE(*,'(A4,3F10.5)') "c3=",m%e(i)%tetC(j,7:9)
            !WRITE(*,'(A4,3F10.5)') "c4=",m%e(i)%tetC(j,10:12)
!
!           Get tet normals -> tet points sequence is changed,
!           so the orientation is the same as in VTK mesh elements
!       
            CALL CalNormal3p( xtri(4,:) , xtri(2,:) , xtri(1,:) , m%e(i)%tetN(j,7:9) )
            CALL CalNormal3p( xtri(4,:) , xtri(3,:) , xtri(2,:) , m%e(i)%tetN(j,1:3) )
            CALL CalNormal3p( xtri(3,:) , xtri(4,:) , xtri(1,:) , m%e(i)%tetN(j,4:6) )
            CALL CalNormal3p( xtri(3,:) , xtri(1,:) , xtri(2,:) , m%e(i)%tetN(j,10:12) )     
          END DO

          npoly = npoly + 1

        END IF
      END DO

!-------------------------------------------------------------------- 
!      ntri = size(m%e(1)%triCon(:,1))
!      DO j = 1,ntri
!        WRITE(*,'(A1,I5,12F6.3)') "C",j,m%e(i)%tetC(j,:)
!        WRITE(*,'(A1,I5,12F6.3)') "N",j,m%e(i)%tetN(j,:)
!      END DO 
!
!      PRINT *, "stopping in calc cent nor" 
!      CALL stopProgram(1)
!--------------------------------------------------------------------

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE GetPolyWeights(m,elid,v,lambda_i)
!
!     Mean value coordinates method for in polyhedron interpolation
!     m --> TYPE(VTKmeshType)
!     elid --> element ID [m%e(elid)]
!     v --> point in which to interpolate nodal values
!     
! -----------------------------------------------------------------------------------------
      USE mesh
      USE logFile      
      IMPLICIT NONE 
      TYPE(VTKmeshType) :: m
  
      INTEGER elid
      REAL(8) v(3)

      INTEGER j,k,l,h
      INTEGER imax, iSide
      
      INTEGER ntri,nv
      INTEGER c1,c2,c3
      REAL(8) xtri(3,3),ui(3),r(3),ri

      REAL(8) wi(m%e(elid)%nVertex)
      REAL(8) lambda_i(m%e(elid)%nVertex)

      REAL(8) rr(m%e(elid)%nVertex,3)

      REAL(8), ALLOCATABLE :: uiMatrix(:,:)
      
      ntri=size(m%e(elid)%triCon(:,1))
      nv=m%e(elid)%nVertex

      ALLOCATE(uiMatrix(nv,ntri+1))

      uiMatrix = 0.0D0
      uiMatrix(:,1) = m%e(elid)%polyUniver

      !DO j = 1,size(uiMatrix(:,1))
      !  PRINT *, uiMatrix(j,1)
      !END DO

      !ntri = 18
      DO j=1,ntri

        IF ( size(m%e(elid)%triCon(1,:)) .EQ. 4 ) THEN

          c1 = m%e(elid)%triCon(j,1)
          c2 = m%e(elid)%triCon(j,2)
          c3 = m%e(elid)%triCon(j,3)

          xtri(1,:) = m%x(c1,:)
          xtri(2,:) = m%x(c2,:)
          xtri(3,:) = m%x(c3,:)

        ELSE IF ( size(m%e(elid)%triCon(1,:)) .EQ. 5 ) THEN

          c1 = m%e(elid)%triCon(j,1)
          c2 = m%e(elid)%triCon(j,2)
          c3 = -2

          xtri(1,:) = m%x(c1,:)
          xtri(2,:) = m%x(c2,:)

          !face center of face "iSide"
          iSide = m%e(elid)%triCon(j,5)
          xtri(3,:) = m%e(elid)%centers(iSide,:)
          
        END IF

        !PRINT *, "writing out triangle weights for j =",j
        CALL GetTriangleWeights(xtri,v,ui)        
        !WRITE(*,'(A4,3F10.5)') "  v=",v
        !WRITE(*,'(A4,I8,F10.5)') "c1=",c1,ui(1)
        !WRITE(*,'(A4,I8,F10.5)') "c2=",c2,ui(2)
        !WRITE(*,'(A4,I8,F10.5)') "c3=",c3,ui(3)

        DO k=1,nv

          IF (c1.EQ.uiMatrix(k,1)) THEN
            uiMatrix(k,j+1) = ui(1)
          ELSE IF (c2.EQ.uiMatrix(k,1)) THEN
            uiMatrix(k,j+1) = ui(2)
          ELSE IF (c3.EQ.uiMatrix(k,1)) THEN
            uiMatrix(k,j+1) = ui(3)
          END IF

        END DO

      END DO      

      DO k=1,nv
        !WRITE(*,'(18F5.2)') uiMatrix(k,:)
        r = m%x(int(uiMatrix(k,1)),:)-v
        CALL vecLen(r,ri)
        wi(k) = sum(uiMatrix(k,2:ntri+1))/ri
      END DO 

      !PRINT *, "Writing out lambda_i"
      DO k=1,nv
        lambda_i(k) = wi(k)/sum(wi)
        !WRITE(*,'(F12.5)') lambda_i(k)
      END DO

      !DO k=1,nv
        !WRITE(*,'(2F12.5)') lambda_i(k), un(m%e(elid)%polyUniver(k))
        !rr(k,1) = lambda_i(k)*un(m%e(elid)%polyUniver(k))
      !END DO
      !uv = sum(rr(:,1))

      !--------------- DEBUG OUTPUT ---------------
      IF (.FALSE.) THEN
        DO k=1,nv
          WRITE(*,'(F9.2)') lambda_i(k)
        END DO

        PRINT *, "test sum: lambda_i(vi-v) -> must be 0"
        DO k=1,nv
          rr(k,:) = lambda_i(k)*(m%x(m%e(elid)%polyUniver(k),:)-v)
        END DO
        WRITE(*,'(3F12.5)') sum(rr(:,1)),sum(rr(:,2)),sum(rr(:,3))

        PRINT *, "test sum: lambda_i*vi -> must be v"
        DO k=1,nv
          WRITE(*,'(4F12.5)') lambda_i(k), m%x(m%e(elid)%polyUniver(k),:)
          rr(k,:) = lambda_i(k)*m%x(m%e(elid)%polyUniver(k),:)
        END DO

        !PRINT *, "test sum: lambda_i*f(vi) -> must be f(v)"
        !DO k=1,nv
        !  WRITE(*,'(5F10.4)') lambda_i(k), un(m%e(elid)%polyUniver(k)), m%x(m%e(elid)%polyUniver(k),:)
        !  rr(k,:) = lambda_i(k)*m%x(m%e(elid)%polyUniver(k),:)
        !END DO

        PRINT *, "-----"
        WRITE(*,'(3F12.5)') sum(rr(:,1)),sum(rr(:,2)),sum(rr(:,3))
        PRINT *, "-----"
      END IF
      !--------------- DEBUG OUTPUT --------------- 

      DEALLOCATE(uiMatrix)

      END 

! -----------------------------------------------------------------------------------------
      SUBROUTINE PolyVolCenterVals(m,fluid)
!
!     Get values at polyhedron volumetric center
!      
! -----------------------------------------------------------------------------------------
      
      USE mesh
      USE superE
      USE mFluid
     
      IMPLICIT NONE 

      TYPE(VTKmeshType) :: m  ! mesh data structure
      TYPE(fluidType) :: fluid ! fluid flow fields

      INTEGER i,j,node

      DO i=1,m%nelem

        IF(m%e(i)%type.EQ.42) THEN

          m%e(i)%vcU(:)=0.0D0
          m%e(i)%vcVort(:)=0.0D0
          m%e(i)%vcGradU(:,:)=0.0D0
          m%e(i)%vcdVdT(:)=0.0D0

          DO j=1,m%e(i)%nVertex

            node = m%e(i)%polyUniver(j)

            !velocity at center
            m%e(i)%vcU = m%e(i)%vcU + fluid%Un(node,:) 

            !vorticity at center
            m%e(i)%vcVort = m%e(i)%vcVort + fluid%Vortn(node,:)

            !velocity gradient at center
            m%e(i)%vcGradU(1,:) =  m%e(i)%vcGradU(1,:) + fluid%gradUxn(node,:)
            m%e(i)%vcGradU(2,:) =  m%e(i)%vcGradU(2,:) + fluid%gradUyn(node,:)
            m%e(i)%vcGradU(3,:) =  m%e(i)%vcGradU(3,:) + fluid%gradUzn(node,:)

            !dV/dT at center
            m%e(i)%vcdVdT(1) = m%e(i)%vcdVdT(1) + fluid%dvxdt(node)
            m%e(i)%vcdVdT(2) = m%e(i)%vcdVdT(2) + fluid%dvydt(node)
            m%e(i)%vcdVdT(3) = m%e(i)%vcdVdT(3) + fluid%dvzdt(node)

          END DO

          m%e(i)%vcU = m%e(i)%vcU / m%e(i)%nVertex
          m%e(i)%vcVort = m%e(i)%vcVort / m%e(i)%nVertex
          m%e(i)%vcGradU = m%e(i)%vcGradU / m%e(i)%nVertex
          m%e(i)%vcdVdT = m%e(i)%vcdVdT / m%e(i)%nVertex

        END IF

      END DO

      END



! -----------------------------------------------------------------------------------------
      SUBROUTINE triangulatePolySurface(m,mode)
!
!
!     TRIANGULATION MODE 1 m%e(i)%triCon(nTri,4) - v1,v2,-1,-1
!     TRIANGULATION MODE 2 m%e(i)%triCon(nTri,5) - v1,v2,-1,-1, face number
!
! -----------------------------------------------------------------------------------------
      USE mesh
      !USE logFile      
      IMPLICIT NONE 

      TYPE(VTKmeshType) :: m

      INTEGER mode
      INTEGER i,j,k,k2,l,l2,ll
      INTEGER ibest,np,ns

      INTEGER nTri,nTriOld
      INTEGER itri_start, itri_end

      REAL(8), ALLOCATABLE :: pnts(:,:)
      REAL(8), ALLOCATABLE :: Fcents(:,:)

      INTEGER, ALLOCATABLE :: conn(:),tricon(:,:)


      !initalize counters for tet decom routine stat
      m%iff=0
      m%inf=0
      m%imf=0

      DO i = 1,m%nelem
        
        IF (m%e(i)%type.EQ.42) THEN

!-----------------------------------------------------
!           TRIANGULATION MODE 1  
!
!           create triangles using only face vertices,
!           using max(minimum angle) principle 
!-----------------------------------------------------
          IF (mode .EQ. 1) THEN

            !-count number of triangles for polyhedron
            nTri = 0
            nTriOld = 0
            l=1
            ns = size(m%e(i)%con(:,1))
            !PRINT *, "Ns=", ns
            DO j=1,ns
              l=l+1
              np = m%e(i)%con(j,1)
              !PRINT *, "np=", np
              nTri = (np -2)+nTriOld
              nTriOld = nTri
              l=l+np
            END DO

            !print *, "i=",i,"ntri=",nTri
            !-end count number of triangles for polyhedron

            ALLOCATE(m%e(i)%triCon(nTri,4))

            itri_start = 0
            itri_end = 0
            DO j=1,ns
              np = m%e(i)%con(j,1)

              ALLOCATE(conn(np),pnts(np,3),tricon(np-2,3))

              !PRINT *, "side=",j

              DO k=1,np
                conn(k) = m%e(i)%con(j,k+1)
                pnts(k,:) = m%x(conn(k),:)
              END DO
              
              !PRINT *, "conn=", conn

              IF (np .GT. 3) THEN

                CALL calcAngles(pnts,np,ibest)
                CALL makeTriangles(pnts,conn,np,ibest,tricon)

              ELSE 

                tricon(1,:) = conn

              END IF

              itri_start = itri_end+1
              itri_end = itri_start+(np-3)

              !PRINT *, itri_start, itri_end

              m%e(i)%triCon(itri_start:itri_end,1:3) = tricon
              m%e(i)%triCon(itri_start:itri_end,4) = -1 !poly element center
              
              DEALLOCATE(conn,pnts,tricon)		

            END DO

!-----------------------------------------------------
!           TRIANGULATION MODE 2 
!
!           triangulation using addition vertex
!           on face centre
!-----------------------------------------------------
          ELSE IF (mode .EQ. 2) THEN

            nTri = 0

            ns=size( m%e(i)%con(:,1) )
            DO j=1,ns
              !WRITE(*,'(7I4)') m%e(i)%con(j,:)
              nTri = nTri + m%e(i)%con(j,1)
            END DO
            !PRINT *, "after smth"

            ALLOCATE(m%e(i)%triCon(nTri,5))

            m%e(i)%triCon(:,:) = -1

            !face center points
            !m%e(i)%centers(j,:)

            l = 0
            DO j=1,ns

              np=m%e(i)%con(j,1)

              IF (np .GT. 3) THEN

                DO k = 1,np
                  l=l+1

                  m%e(i)%triCon(l,1) = m%e(i)%con(j,k+1)

                  IF (k .LT. np) THEN
                    m%e(i)%triCon(l,2) = m%e(i)%con(j,k+2)
                  ELSE
                    m%e(i)%triCon(l,2) = m%e(i)%con(j,2)
                  END IF

                  m%e(i)%triCon(l,3) = -1
                  m%e(i)%triCon(l,4) = -1

                  !indicate to which side it belongs
                  m%e(i)%triCon(l,5) = j
                END DO

              ELSE 

                l=l+1
                m%e(i)%triCon(l,1:3) = m%e(i)%con(j,2:4)
                m%e(i)%triCon(l,4) = -1

              END IF

            END DO

          END IF

          !write out triCon
          !PRINT *, "writ triCon"
          !DO j=1,size(m%e(i)%triCon(:,1))
          !	WRITE(*,'(5I4)') m%e(i)%triCon(j,:)
          !END DO
          !PRINT *, "end writ triCon"

        END IF

      END DO


      !DO i = 1,m%nelem

        !Print the first polyhedron element
        !IF (m%e(i)%type.EQ.42) THEN

        !Print specific element number, if you want to inspect specific
      !  IF (i.EQ.1720) THEN
      !    PRINT *, "Printing out tricon for element:",i
      !    DO j = 1,size( m%e(i)%triCon(:,1)  )
      !      PRINT *, m%e(i)%triCon(j,:)
      !    END DO
      !    EXIT
      !  END IF
      !END DO

      DO i = 1,m%nelem
        
        !Print the first polyhedron element
        !IF (m%e(i)%type.EQ.42) THEN

        !Print specific element number, if you want to inspect specific
        IF (i.EQ.1720) THEN

          PRINT *, "Printing out faces comparison for element:",i
         
          DO j = 1,size(m%e(i)%con(:,1))
            PRINT *,  "Econ ",j,m%e(i)%con(j, 2:m%e(i)%con(j,1)+1  )
          END DO
          
          EXIT
        
        END IF
      END DO

      PRINT *, "nfac=",m%nfac

      k=0
      DO j = 1,m%nfac
        IF (m%f(j)%owner.EQ.i) THEN
          k=k+1
          !PRINT *, "Ownr=",m%f(j)%owner
          !PRINT *, "Neib=",m%f(j)%neighbour
          PRINT *, "FconO=",k,m%f(j)%con(:)        
        END IF
        IF (m%f(j)%neighbour.EQ.i) THEN
          k=k+1
          !PRINT *, "Ownr=",m%f(j)%owner
          !PRINT *, "Neib=",m%f(j)%neighbour
          PRINT *, "FconN=",k,m%f(j)%con(:)        
        END IF
      END DO




      !PRINT *, "STOPPING IN triang"
      !CALL stopProgram(1)

      END

! -----------------------------------------------------------------------------------------
      SUBROUTINE GetTriangleWeights(x,v,ui)     
! -----------------------------------------------------------------------------------------
    
      IMPLICIT NONE 
      REAL(8) x(3,3), v(3)

      INTEGER i,j,k
      REAL(8) vi(3),vj(3),vk(3)
      REAL(8) ei(3),ej(3),ek(3)

      REAL(8) bij,bjk,bki
      REAL(8) nij(3),njk(3),nki(3)

      REAL(8) dp1,dp2,dp3
      
      REAL(8) ui(3)

      LOGICAL correcting

      correcting = .FALSE.
      
10    IF(correcting) v = v + 1.0E-08   

      DO i = 1,3

        IF (i.EQ.1) THEN
          vi = x(1,:)-v
          vj = x(2,:)-v
          vk = x(3,:)-v
        ELSE IF (i.EQ.2) THEN
          vi = x(2,:)-v
          vj = x(3,:)-v
          vk = x(1,:)-v
        ELSE IF (i.EQ.3) THEN
          vi = x(3,:)-v
          vj = x(1,:)-v
          vk = x(2,:)-v
        END IF

        CALL vecAngle(vi,vj,bij)
        CALL vecAngle(vj,vk,bjk)
        CALL vecAngle(vk,vi,bki)

        CALL CalNormal2vec(vi,vj,nij)
        CALL CalNormal2vec(vj,vk,njk)
        CALL CalNormal2vec(vk,vi,nki)

        ei = vi
        ej = vj
        ek = vk

        CALL Normalize(ei)
        CALL Normalize(ej)
        CALL Normalize(ek)
        
        CALL DotProduct(nij,njk,dp1)
        CALL DotProduct(nki,njk,dp2)
        CALL DotProduct(ei,njk,dp3)

        IF (correcting) THEN
!          WRITE(*,'(A)') "-------------------------------------------"
!          WRITE(*,'(A)') "After:"
!
!          PRINT *, "writing out"
!          WRITE(*,'(A4,3F10.5)') "  v=",v
!
!          WRITE(*,'(A6,3F10.5)') "x(1,:)=",x(1,:)
!          WRITE(*,'(A6,3F10.5)') "x(2,:)=",x(2,:)
!          WRITE(*,'(A6,3F10.5)') "x(3,:)=",x(3,:)
!
!          WRITE(*,'(A4,3F10.5)') " ei=",ei
!          WRITE(*,'(A4,3F10.5)') " ej=",ej
!          WRITE(*,'(A4,3F10.5)') " ek=",ek
!
!          WRITE(*,'(A4,3F10.5)') "njk=",njk
!          WRITE(*,'(A4,3F10.5)') "dp3=",dp3
          correcting = .FALSE.
        END IF

        IF (ABS(dp3).LT.1.0E-08) THEN
!          WRITE(*,'(A4,F10.5,A22)') "dp3=",dp3,"correcting"
!
!          PRINT *, "writing out"
!          WRITE(*,'(A4,3F10.5)') "  v=",v
!
!          WRITE(*,'(A6,3F10.5)') "x(1,:)=",x(1,:)
!          WRITE(*,'(A6,3F10.5)') "x(2,:)=",x(2,:)
!          WRITE(*,'(A6,3F10.5)') "x(3,:)=",x(3,:)
!
!          WRITE(*,'(A4,3F10.5)') " ei=",ei
!          WRITE(*,'(A4,3F10.5)') " ej=",ej
!          WRITE(*,'(A4,3F10.5)') " ek=",ek
!
!          WRITE(*,'(A4,3F10.5)') "njk=",njk
!          WRITE(*,'(A4,3F10.5)') "dp3=",dp3
!          !WRITE(*,'(A4,I8,F10.5)') "c1=",c1,ui(1)
!          !WRITE(*,'(A4,I8,F10.5)') "c2=",c2,ui(2)
!          !WRITE(*,'(A4,I8,F10.5)') "c3=",c3,ui(3)
!
          correcting=.TRUE.
          GOTO 10
        END IF

        ui(i) = (bjk + bij*dp1 + bki*dp2) / (2*dp3)

      END DO

      !DO i = 1,3
      !  PRINT *, "ui=", ui(i)
      !END DO


      END


! -----------------------------------------------------------------------------------------
      SUBROUTINE makeTriangles(pnts,con,n,ibest,tricon)
!
! -----------------------------------------------------------------------------------------

      IMPLICIT NONE 
      INTEGER n,nTri,nCon,i,j,k,l,ibest
      INTEGER con(n),newcon(n),tricon(n-2,3)
      REAL(8) pnts(n,3)

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