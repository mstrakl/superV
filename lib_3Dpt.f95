!
!
!     lib_3Dpt  - 3D particle tracking
!
!     written by Jure Ravnik (jure.ravnik@um.si)
!
!     List of  Routines
!                      - hex_kks2lks  ! find xi,eta,zeta from x,y,z for hex element
!                      - hex_lks2kks  ! find x,y,z from xi,eta,zeta for hex element
!                      - tet_kks2lks  ! find xi,eta,zeta from x,y,z for tet element
!                      - pri_kks2lks  ! find xi,eta,zeta from x,y,z for prizm element
!                      - pyr_kks2lks  ! find xi,eta,zeta from x,y,z for pyramid element
!
!     Version : 14. November 2018
!

! ------------------------------------------------------------------------------------    
      SUBROUTINE tet_kks2lks(r,p1,p2,p3,p4,xez,ierr)
!
!     VTK distribution of nodes, Cramers rule for solution of 3x3 SLE
!
!           4
!          / \  
!         /   \   3
!        /     \ /
!       1 ----- 2
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) xez(3),r(3),p1(3),p2(3),p3(3),p4(3)
      REAL(8) c1(3),c2(3),c3(3),b(3) ! sys matrix colums
      REAL(8) d1,d2,d3,d ! determinants
      INTEGER i,ierr

      DO i=1,3
        c1(i)=p2(i)-p1(i)
        c2(i)=p3(i)-p1(i)
        c3(i)=p4(i)-p1(i)
        b(i) = r(i)-p1(i)
      END DO
      
      CALL CalDeterminant(c1,c2,c3,d)
      IF (d.EQ.0.0D0) THEN
        ierr=1
        RETURN
      END IF

      CALL CalDeterminant(b,c2,c3,d1)
      CALL CalDeterminant(c1,b,c3,d2)
      CALL CalDeterminant(c1,c2,b,d3)

      xez(1)=d1/d
      xez(2)=d2/d
      xez(3)=d3/d
     
      ierr = 0

      END

! ------------------------------------------------------------------------------------ 
      SUBROUTINE SetMatrixInverse3x3(Jac,iJac)
      IMPLICIT NONE
      REAL(8) Jac(3,3),iJac(3,3),ajac

!
!     Jacobi matrix determinant
!
      CALL CalDeterminantMatrix(Jac,ajac)
!
!     Invert Jacobi matrix -> d(xi)/d(x)
!
!
!     prvi stolpec
      iJac(1,1) = (Jac(3,3)*Jac(2,2)-Jac(2,3)*Jac(3,2))/ajac ! dxidx
      iJac(2,1) = (Jac(2,3)*Jac(3,1)-Jac(2,1)*Jac(3,3))/ajac ! dxidy
      iJac(3,1) = (Jac(2,1)*Jac(3,2)-Jac(2,2)*Jac(3,1))/ajac ! dxidz
!     drugi stolpec           
      iJac(1,2) = (Jac(1,3)*Jac(3,2)-Jac(3,3)*Jac(1,2))/ajac ! detadx
      iJac(2,2) = (Jac(3,3)*Jac(1,1)-Jac(1,3)*Jac(3,1))/ajac ! detady
      iJac(3,2) = (Jac(1,2)*Jac(3,1)-Jac(1,1)*Jac(3,2))/ajac ! detadz
!     tretji stolpec           
      iJac(1,3) = (Jac(2,3)*Jac(1,2)-Jac(1,3)*Jac(2,2))/ajac ! dzetadx
      iJac(2,3) = (Jac(2,1)*Jac(1,3)-Jac(2,3)*Jac(1,1))/ajac ! dzetady
      iJac(3,3) = (Jac(1,1)*Jac(2,2)-Jac(2,1)*Jac(1,2))/ajac ! dzetadz

      END


! ------------------------------------------------------------------------------------ 
      SUBROUTINE CalDeterminantMatrix(c,d)
      IMPLICIT NONE
      REAL(8) c(3,3),d

      d = c(1,1)*c(2,2)*c(3,3) + &
     &    c(1,2)*c(2,3)*c(3,1) +       &
     &    c(1,3)*c(2,1)*c(3,2) -            &
     &    c(3,1)*c(2,2)*c(1,3) -            &
     &    c(1,1)*c(3,2)*c(2,3) -            &
         c(2,1)*c(1,2)*c(3,3)    

      END


! ------------------------------------------------------------------------------------ 
      SUBROUTINE CalDeterminant(c1,c2,c3,d)
      IMPLICIT NONE
      REAL(8) c1(3),c2(3),c3(3),d

      d = c1(1)*c2(2)*c3(3) + &
     &    c2(1)*c3(2)*c1(3) +       &
     &    c3(1)*c1(2)*c2(3) -            &
     &    c1(3)*c2(2)*c3(1) -            &
     &    c1(1)*c2(3)*c3(2) -            &
         c1(2)*c2(1)*c3(3)           

      END

! ------------------------------------------------------------------------------------    
      SUBROUTINE tet_shapef(xez,fig)
!
!     VTK distribution of nodes
!
!           4
!          / \  
!         /   \   3
!        /     \ /
!       1 ----- 2
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) xez(3),fig(4)

      FIG(1)=1.0D0-xez(1)-xez(2)-xez(3)
      FIG(2)=xez(1)
      FIG(3)=xez(2)
      FIG(4)=xez(3)


      END

! ------------------------------------------------------------------------------------    
      SUBROUTINE tet_shapef_der(xez,SFD)

      IMPLICIT NONE
      REAL(8) xez(3),SFD(4,3)

      SFD(1,1) = -1.0D0       ! d( FIG(1) )/d(xi)
      SFD(1,2) = -1.0D0       ! d( FIG(1) )/d(eta)
      SFD(1,3) = -1.0D0       ! d( FIG(1) )/d(zeta)

      SFD(2,1) = +1.0D0       ! d( FIG(4) )/d(xi)
      SFD(2,2) =  0.0D0       ! d( FIG(4) )/d(eta)
      SFD(2,3) =  0.0D0      ! d( FIG(4) )/d(zeta)

      SFD(3,1) =  0.0D0          ! d( FIG(3) )/d(xi)
      SFD(3,2) = +1.0D0        ! d( FIG(3) )/d(eta)
      SFD(3,3) =  0.0D0      ! d( FIG(3) )/d(zeta)

      SFD(4,1) =  0.0D0         ! d( FIG(2) )/d(xi)
      SFD(4,2) =  0.0D0          ! d( FIG(2) )/d(eta)
      SFD(4,3) = +1.0D0     ! d( FIG(2) )/d(zeta)

      END

! ------------------------------------------------------------------------------------    
      SUBROUTINE pyr_shapef(xez,FIG)
!
!     VTK distribution of nodes
!
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
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) xez(3),FIG(5),fac

      IF (xez(3).EQ.1.0D0) THEN
        fac=-xez(3)
      ELSE
        fac=-xez(3)+xez(1)*xez(2)*xez(3)/(1.0D0-xez(3))
      END IF

      FIG(1)=0.25D0*( (1.0D0-xez(1))*(1.0D0-xez(2)) + fac )
      FIG(2)=0.25D0*( (1.0D0+xez(1))*(1.0D0-xez(2)) + fac )
      FIG(3)=0.25D0*( (1.0D0+xez(1))*(1.0D0+xez(2)) + fac )
      FIG(4)=0.25D0*( (1.0D0-xez(1))*(1.0D0+xez(2)) + fac )
      FIG(5)=xez(3)

      END

! ------------------------------------------------------------------------------------    
      SUBROUTINE pyr_shapef_der(xez,SFD)

      IMPLICIT NONE
      REAL(8) xez(3),SFD(5,3),fac

      IF (xez(3).EQ.1.0D0) THEN
        fac=0.0D0
      ELSE
        fac=1.0D0/(1.0D0-xez(3))
      END IF

      SFD(1,1) = 0.25D0 * ( -1.0D0 + xez(2) + xez(3)*xez(2)*fac )      ! d( FIG(1) )/d(xi)
      SFD(1,2) = 0.25D0 * ( -1.0D0 + xez(1) + xez(3)*xez(1)*fac )      ! d( FIG(1) )/d(eta)
      SFD(1,3) = 0.25D0 * ( -1.0D0 + xez(2)*xez(1)*fac + xez(1)*xez(2)*xez(3)*fac*fac )  ! d( FIG(1) )/d(zeta)

      SFD(2,1) = 0.25D0 * ( +1.0D0 - xez(2) + xez(3)*xez(2)*fac )      ! d( FIG(4) )/d(xi)
      SFD(2,2) = 0.25D0 * ( -1.0D0 - xez(1) + xez(3)*xez(1)*fac )      ! d( FIG(4) )/d(eta)
      SFD(2,3) = 0.25D0 * ( -1.0D0 + xez(2)*xez(1)*fac + xez(1)*xez(2)*xez(3)*fac*fac )   ! d( FIG(4) )/d(zeta)

      SFD(3,1) = 0.25D0 * ( +1.0D0 + xez(2) + xez(3)*xez(2)*fac )      ! d( FIG(3) )/d(xi)
      SFD(3,2) = 0.25D0 * ( +1.0D0 + xez(1) + xez(3)*xez(1)*fac )      ! d( FIG(3) )/d(eta)
      SFD(3,3) = 0.25D0 * ( -1.0D0 + xez(2)*xez(1)*fac + xez(1)*xez(2)*xez(3)*fac*fac )  ! d( FIG(3) )/d(zeta)

      SFD(4,1) = 0.25D0 * ( -1.0D0 - xez(2) + xez(3)*xez(2)*fac )      ! d( FIG(2) )/d(xi)
      SFD(4,2) = 0.25D0 * ( +1.0D0 - xez(1) + xez(3)*xez(1)*fac )      ! d( FIG(2) )/d(eta)
      SFD(4,3) = 0.25D0 * ( -1.0D0 + xez(2)*xez(1)*fac + xez(1)*xez(2)*xez(3)*fac*fac )  ! d( FIG(2) )/d(zeta)

      SFD(5,1) = 0.00D0         ! d( FIG(5) )/d(xi)
      SFD(5,2) = 0.00D0         ! d( FIG(5) )/d(eta)
      SFD(5,3) = 1.00D0         ! d( FIG(5) )/d(zeta)


      END


! ------------------------------------------------------------------------------------    
      SUBROUTINE pyr_FindZeta(r,n,c,p,zeta)
!
!     (ne dela v splosnem ) 
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3) ! point inside element
      REAL(8) p(3) ! fifth point of pyramid
      REAL(8) n(3) ! normal to square face
      REAL(8) c(3) ! center of square face
      REAL(8) zeta ! zeta
      REAL(8) tot

      CALL DotProduct(n,r-c,zeta)
      CALL DotProduct(n,p-c,tot)

      zeta=zeta/tot

      END


! ------------------------------------------------------------------------------------    
      SUBROUTINE pyr_kks2lks(r,p,xez,ierr)

!
!     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
!       minimiziramo in odvode jacobijeve 
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3),rr(3) ! point inside element
      REAL(8) p(5,3) ! vertex locations
      REAL(8) xez(3) ! xi,eta,zeta

      REAL(8) eps
      INTEGER maxit,nit,ierr

!     zacetni priblizek
      xez(1)=0.0D0
      xez(2)=0.0D0
      xez(3)=0.0D0      

!     toleranca 
      eps=1.0D-10 !-6 !-13
!     najvecje stevilo iteracij
      maxit=1000 !10000

      CALL pyr_mnewt(maxit,xez,3,eps,eps,r,p,nit)

!     zapomnimo rezultate, preverimo ce skonvergiralo !!
      ierr=0
      IF (nit.GE.maxit) THEN
!        WRITE (*,*) "Dosezeno maximalno stevilo iteracij, kks2lks",nit
!
!       Use DOE to find xi,eta,zeta since Newton-Raphsen failed
!
        CALL pyr_kks2lks_doe(r,p,xez,100)
        ierr=0
      END IF

      END



! ------------------------------------------------------------------------------------    
      SUBROUTINE pyr_kks2lks_doe(r,p,xez,maxit)
!
!     $: DOE for xi,eta,zeta
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3) ! point inside element
      REAL(8) p(5,3) ! vertex locations
      REAL(8) xez(3),x(3) ! xi,eta,zeta

      REAL(8) dx,rr(3),minD,d
      INTEGER maxit,i,j,k
 
      dx = 2.0D0 / maxit

      minD=1.0D10
      DO i=0,maxit
        x(1)=-1.0D0+dx*i
        DO j=0,maxit
          x(2)=-1.0D0+dx*j
          DO k=0,maxit
            x(3)=0.5D0*dx*k
            CALL pyr_lks2kkk(x,p,rr)
            d=(r(1)-rr(1))**2+(r(2)-rr(2))**2+(r(3)-rr(3))**2
            IF (d.LT.minD) THEN  
              minD=d
              xez=x
            END IF
          END DO
        END DO
      END DO


      END

! ------------------------------------------------------------------------------------    
      SUBROUTINE pyr_lks2kkk(xez,p,r)

!
!     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
!       minimiziramo in odvode jacobijeve 
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3) ! point inside element
      REAL(8) p(5,3) ! vertex locations
      REAL(8) xez(3) ! xi,eta,zeta
      REAL(8) FIG(5) ! shape functions
      INTEGER i,j,k
!
!     Calculate shape functions
!
      CALL pyr_shapef(xez,FIG)
!
!     Functions to be minimized
!
      r=0.0D0
      DO i=1,3 ! x,y,z
        DO j=1,5 ! nodes
          r(i)=r(i)+p(j,i)*FIG(j)
        END DO
      END DO

      END
     

! ------------------------------------------------------------------------------------    
      SUBROUTINE pyr_kks2lks2(r,p,xez,ierr)

!
!     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
!       minimiziramo in odvode jacobijeve 
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3) ! point inside element
      REAL(8) p(5,3) ! vertex locations
      REAL(8) xez(3) ! xi,eta,zeta
      REAL(8) xx(2)

      REAL(8) eps
      INTEGER maxit,nit,ierr

!     zacetni priblizek
      xx(1)=0.0D0
      xx(2)=0.0D0

!     toleranca 
      eps=1.0D-10 !-6 !-13
!     najvecje stevilo iteracij
      maxit=1000 !10000
            
      CALL pyr_mnewt2(maxit,xx,2,eps,eps,r,p,nit,xez(3))

      xez(1)=xx(1)
      xez(2)=xx(2)

!     zapomnimo rezultate, preverimo ce skonvergiralo !!
      ierr=0
      IF (nit.GE.maxit) THEN
!        WRITE (*,*) "Dosezeno maximalno stevilo iteracij, pyr_kks2lks2",nit
        xez(1)=0.0D0
        xez(2)=0.0D0
        ierr=1 
      END IF

      END
               

! ------------------------------------------------------------------------------------    
      SUBROUTINE pyr_kkk2lks_userfun(r,p,xez,F,Jac)

!
!     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
!       minimiziramo in odvode jacobijeve 
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3) ! point inside element
      REAL(8) p(5,3) ! vertex locations
      REAL(8) xez(3) ! xi,eta,zeta
      REAL(8) FIG(5) ! shape functions
      
      REAL(8) F(3),Jac(3,3),SFD(5,3)
      INTEGER i,j,k
!
!     Calculate shape functions
!
      CALL pyr_shapef(xez,FIG)
!
!     Functions to be minimized
!
      DO i=1,3 ! x,y,z
        F(i)=-r(i)
        DO j=1,5 ! nodes
          F(i)=F(i)+p(j,i)*FIG(j)
        END DO
      END DO
!
!     Calculate shape function derivatives
!
      CALL pyr_shapef_der(xez,SFD)
!
!     Jacobi matrix of derivatives J_ij=\p F_i / \p xez_j   
!
      DO i=1,3 ! Fx,Fy,Fz
        DO j=1,3 ! xi,eta,zeta
          Jac(i,j)=0.0D0
          DO k=1,5 ! nodes
            Jac(i,j) = Jac(i,j) + SFD(k,j) * p(k,i)
          END DO
        END DO 
      END DO
       
      END

      
! ------------------------------------------------------------------------------------        
      SUBROUTINE pyr_mnewt(ntrial,xx,n,tolx,tolf,r,points,nit)
      INTEGER n,ntrial
      REAL(8) tolf,tolx,xx(n)
      REAL(8) r(3),points(5,3)

!U    USES lubksb,ludcmp,usrfun
      INTEGER i,k,indx(n),nit
      REAL(8) d,errf,errx,fjac(n,n),fvec(n),p(n)
      do 14  k=1,ntrial
        nit=k
        CALL pyr_kkk2lks_userfun(r,points,xx,Fvec,Fjac)
        errf=0.
        do 11 i=1,n
          errf=errf+abs(fvec(i))
11      continue
        if(errf.le.tolf)return
        do 12 i=1,n
          p(i)=-fvec(i)
12      continue
        call ludcmp(fjac,n,n,indx,d)
        call lubksb(fjac,n,n,indx,p)
        errx=0.
        do 13 i=1,n
          errx=errx+abs(p(i))
          xx(i)=xx(i)+p(i)
13      continue
        if(errx.le.tolx)return
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.  

      
! ------------------------------------------------------------------------------------        
      SUBROUTINE pyr_mnewt2(ntrial,xx,n,tolx,tolf,r,points,nit,zeta)
      INTEGER n,ntrial
      REAL(8) tolf,tolx,xx(n)
      REAL(8) r(3),points(5,3),zeta

!U    USES lubksb,ludcmp,usrfun
      INTEGER i,k,indx(n),nit
      REAL(8) d,errf,errx,fjac(n,n),fvec(n),p(n)
      do 14  k=1,ntrial
        nit=k
        CALL pyr_kkk2lks_userfun2(r,points,xx,Fvec,Fjac,zeta)
        errf=0.
        do 11 i=1,n
          errf=errf+abs(fvec(i))
11      continue
        if(errf.le.tolf)return
        do 12 i=1,n
          p(i)=-fvec(i)
12      continue
        call ludcmp(fjac,n,n,indx,d)
        call lubksb(fjac,n,n,indx,p)
        errx=0.
        do 13 i=1,n
          errx=errx+abs(p(i))
          xx(i)=xx(i)+p(i)
13      continue
        if(errx.le.tolx)return
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.  


! ------------------------------------------------------------------------------------    
      SUBROUTINE pyr_kkk2lks_userfun2(r,p,xx,F,Jac,zeta)

!
!     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
!       minimiziramo in odvode jacobijeve 
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3) ! point inside element
      REAL(8) p(5,3) ! vertex locations
      REAL(8) xx(2) ! xi,eta,zeta
      REAL(8) FIG(5) ! shape functions
      
      REAL(8) F(2),Jac(2,2),SFD(5,3),zeta,xez(3)
      INTEGER i,j,k

      xez(1)=xx(1)
      xez(2)=xx(2)
      xez(3)=zeta   
!
!     Calculate shape functions
!
      CALL pyr_shapef(xez,FIG)
!
!     Functions to be minimized
!
      DO i=1,2 ! x,y,z
        F(i)=-r(i)
        DO j=1,5 ! nodes
          F(i)=F(i)+p(j,i)*FIG(j)
        END DO
      END DO
!
!     Calculate shape function derivatives
!
      CALL pyr_shapef_der(xez,SFD)
!
!     Jacobi matrix of derivatives J_ij=\p F_i / \p xez_j   
!
      DO i=1,2 ! Fx,Fy,Fz
        DO j=1,2 ! xi,eta,zeta
          Jac(i,j)=0.0D0
          DO k=1,5 ! nodes
            Jac(i,j) = Jac(i,j) + SFD(k,j) * p(k,i)
          END DO
        END DO 
      END DO
       
      END




! ------------------------------------------------------------------------------------    
      SUBROUTINE pri_shapef(xez,FIG)
!
!     VTK distribution of nodes
!     0 < xi < 1-eta
!     0 < eta
!     -1 < zeta < 1
!
!           
!                 6         1 (0,0,-1)
!                /|         2 (1,0,-1)
!       4 ----- 5 |         3 (0,1,-1)
!       |       | |         4 (0,0,+1)
!       |       | 3         5 (1,0,+1)
!       |       |/          6 (0,1,+1)
!       1 ----- 2
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) xez(3),FIG(6)

      FIG(1)=0.5D0*(1.0D0-xez(1)-xez(2))*(1.0D0-xez(3))
      FIG(2)=0.5D0*xez(1)*(1.0D0-xez(3))
      FIG(3)=0.5D0*xez(2)*(1.0D0-xez(3))
      FIG(4)=0.5D0*(1.0D0-xez(1)-xez(2))*(1.0D0+xez(3))
      FIG(5)=0.5D0*xez(1)*(1.0D0+xez(3))
      FIG(6)=0.5D0*xez(2)*(1.0D0+xez(3))

      END
! ------------------------------------------------------------------------------------    
      SUBROUTINE pri_shapef_der(xez,SFD)

      IMPLICIT NONE
      REAL(8) xez(3),SFD(6,3)

      SFD(1,1) =   0.5D0 * ( -1.0D0 + xez(3) )          ! d( FIG(1) )/d(xi)
      SFD(1,2) =   0.5D0 * ( -1.0D0 + xez(3) )          ! d( FIG(1) )/d(eta)
      SFD(1,3) =   0.5D0 * ( -1.0D0 + xez(1) + xez(2) ) ! d( FIG(1) )/d(zeta)

      SFD(2,1) =   0.5D0 * ( +1.0D0 - xez(3) )          ! d( FIG(2) )/d(xi)
      SFD(2,2) =   0.0D0                                ! d( FIG(2) )/d(eta)
      SFD(2,3) = - 0.5D0 * xez(1)                       ! d( FIG(2) )/d(zeta)

      SFD(3,1) =   0.0D0                                ! d( FIG(3) )/d(xi)
      SFD(3,2) =   0.5D0 * ( +1.0D0 - xez(3) )          ! d( FIG(3) )/d(eta)
      SFD(3,3) = - 0.5D0 * xez(2)                       ! d( FIG(3) )/d(zeta)

      SFD(4,1) =   0.5D0 * ( -1.0D0 - xez(3) )          ! d( FIG(4) )/d(xi)
      SFD(4,2) =   0.5D0 * ( -1.0D0 - xez(3) )          ! d( FIG(4) )/d(eta)
      SFD(4,3) =   0.5D0 * ( +1.0D0 - xez(1) - xez(2) ) ! d( FIG(4) )/d(zeta)

      SFD(5,1) =   0.5D0 * ( +1.0D0 + xez(3) )          ! d( FIG(5) )/d(xi)
      SFD(5,2) =   0.0D0                                ! d( FIG(5) )/d(eta)
      SFD(5,3) =   0.5D0 * xez(1)                       ! d( FIG(5) )/d(zeta)

      SFD(6,1) =   0.0D0                                ! d( FIG(6) )/d(xi)
      SFD(6,2) =   0.5D0 * ( +1.0D0 + xez(3) )          ! d( FIG(6) )/d(eta)
      SFD(6,3) =   0.5D0 * xez(2)                       ! d( FIG(6) )/d(zeta)

      END


! ------------------------------------------------------------------------------------    
      SUBROUTINE pri_kks2lks_org(r,p,xez,ierr)

!
!     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
!       minimiziramo in odvode jacobijeve 
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3) ! point inside element
      REAL(8) p(6,3) ! vertex locations
      REAL(8) xez(3) ! xi,eta,zeta

      REAL(8) eps
      INTEGER maxit,nit,ierr

!     zacetni priblizek
      xez=0.0D0
!     toleranca 
      eps=1.0D-10 !-6 !-13
!     najvecje stevilo iteracij
      maxit=1000 !10000
            
      CALL pri_mnewt(maxit,xez,3,eps,eps,r,p,nit)

!     zapomnimo rezultate, preverimo ce skonvergiralo !!
      IF (nit.GE.maxit) THEN
!        WRITE (*,*) "Dosezeno maximalno stevilo iteracij, pri_kks2lks",nit
        ierr=1 
        RETURN
      END IF

!     preverimo, ce je noter ali zunaj
      IF (xez(1).GE.-1.0D0.AND.xez(1).LE.1.0D0.AND. &
     &    xez(2).GE.-1.0D0.AND.xez(2).LE.1.0D0.AND. &
         xez(3).GE.-1.0D0.AND.xez(3).LE.1.0D0) THEN
        ierr=0
      ELSE 
        ierr=1 
      END IF

      END

! ------------------------------------------------------------------------------------    
      SUBROUTINE pri_kks2lks(r,p,xez,ierr)

!
!     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
!       minimiziramo in odvode jacobijeve 
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3) ! point inside element
      REAL(8) p(6,3) ! vertex locations
      REAL(8) xez(3) ! xi,eta,zeta

      REAL(8) eps,minVal,maxVal
      INTEGER maxit,nit,ierr

!     zacetni priblizek
      xez=0.0D0
!     toleranca 
      eps=1.0D-10 !-6 !-13
!     najvecje stevilo iteracij
      maxit=1000 !10000
            
      CALL pri_mnewt(maxit,xez,3,eps,eps,r,p,nit)

!     zapomnimo rezultate, preverimo ce skonvergiralo !!
      IF (nit.GE.maxit) THEN
!        WRITE (*,*) "Dosezeno maximalno stevilo iteracij, pri_kks2lks",nit
        ierr=1 
        RETURN
      END IF

      minVal = -1.0D0 - 5e-01
      maxVal = 1.0D0 + 5e-01

!     preverimo, ce je noter ali zunaj
      IF (xez(1).GE.minVal.AND.xez(1).LE.maxVal.AND. &
     &    xez(2).GE.minVal.AND.xez(2).LE.maxVal.AND. &
         xez(3).GE.minVal.AND.xez(3).LE.maxVal) THEN
        ierr=0
      ELSE 
        ierr=1 
      END IF


      END      

! ------------------------------------------------------------------------------------    
      SUBROUTINE pri_kkk2lks_userfun(r,p,xez,F,Jac)

!
!     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
!       minimiziramo in odvode jacobijeve 
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) r(3) ! point inside element
      REAL(8) p(6,3) ! vertex locations
      REAL(8) xez(3) ! xi,eta,zeta
      REAL(8) FIG(6) ! shape functions
      
      REAL(8) F(3),Jac(3,3),SFD(6,3)
      INTEGER i,j,k
!
!     Calculate shape functions
!
      CALL pri_shapef(xez,FIG)
!
!     Functions to be minimized
!
      DO i=1,3 ! x,y,z
        F(i)=-r(i)
        DO j=1,6 ! nodes
          F(i)=F(i)+p(j,i)*FIG(j)
        END DO
      END DO
!
!     Calculate shape function derivatives
!
      CALL pri_shapef_der(xez,SFD)
!
!     Jacobi matrix of derivatives J_ij=\p F_i / \p xez_j   
!
      DO i=1,3 ! Fx,Fy,Fz
        DO j=1,3 ! xi,eta,zeta
          Jac(i,j)=0.0D0
          DO k=1,6 ! nodes
            Jac(i,j) = Jac(i,j) + SFD(k,j) * p(k,i)
          END DO
        END DO 
      END DO
       
      END

      
! ------------------------------------------------------------------------------------        
      SUBROUTINE pri_mnewt(ntrial,xx,n,tolx,tolf,r,points,nit)
      INTEGER n,ntrial
      REAL(8) tolf,tolx,xx(n)
      REAL(8) r(3),points(6,3)

!U    USES lubksb,ludcmp,usrfun
      INTEGER i,k,indx(n),nit
      REAL(8) d,errf,errx,fjac(n,n),fvec(n),p(n)
      do 14  k=1,ntrial
        nit=k
        CALL pri_kkk2lks_userfun(r,points,xx,Fvec,Fjac)
        errf=0.
        do 11 i=1,n
          errf=errf+abs(fvec(i))
11      continue
        if(errf.le.tolf)return
        do 12 i=1,n
          p(i)=-fvec(i)
12      continue
        call ludcmp(fjac,n,n,indx,d)
        call lubksb(fjac,n,n,indx,p)
        errx=0.
        do 13 i=1,n
          errx=errx+abs(p(i))
          xx(i)=xx(i)+p(i)
13      continue
        if(errx.le.tolx)return
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.  


! ------------------------------------------------------------------------------------    
      SUBROUTINE hex_shapef(xez,fig)
!
!     VTK distribution of nodes
!
!
!          8 ---- 7
!         /|     /|
!        / |    / |
!       5 ---- 6  |
!       |  |   |  |
!       |  4 --|- 3
!       | /    |/
!       1 ---- 2
!      
! ------------------------------------------------------------------------------------    
      IMPLICIT NONE
      REAL(8) xez(3),fig(8),eta1m,eta1p,eta2m,eta2p,eta3m,eta3p

      ETA1M=1.0D0-xez(1)
      ETA1P=1.0D0+xez(1)
      ETA2M=1.0D0-xez(2)
      ETA2P=1.0D0+xez(2) 
      ETA3M=1.0D0-xez(3)  
      ETA3P=1.0D0+xez(3) 

      FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
      FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
      FIG(4)=0.125D0*ETA1M*ETA2P*ETA3M
      FIG(3)=0.125D0*ETA1P*ETA2P*ETA3M
      FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
      FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
      FIG(8)=0.125D0*ETA1M*ETA2P*ETA3P
      FIG(7)=0.125D0*ETA1P*ETA2P*ETA3P

      END


! ------------------------------------------------------------------------------------    
      SUBROUTINE hex_shapef_der(xez,SFD)

      IMPLICIT NONE
      REAL(8) xez(3),SFD(8,3)

!
!          8 ---- 7
!         /|     /|
!        / |    / |
!       5 ---- 6  |
!       |  |   |  |
!       |  4 --|- 3
!       | /    |/
!       1 ---- 2

      SFD(1,1) = - 0.125D0 * (1.0D0 - xez(3)) * (1.0D0 - xez(2))     ! d( FIG(1) )/d(xi)
      SFD(1,2) = - 0.125D0 * (1.0D0 - xez(3)) * (1.0D0 - xez(1))     ! d( FIG(1) )/d(eta)
      SFD(1,3) = - 0.125D0 * (1.0D0 - xez(2)) * (1.0D0 - xez(1))     ! d( FIG(1) )/d(zeta)

      SFD(2,1) =   0.125D0 * (1.0D0 - xez(3)) * (1.0D0 - xez(2))     ! d( FIG(2) )/d(xi)
      SFD(2,2) = - 0.125D0 * (1.0D0 - xez(3)) * (1.0D0 + xez(1))     ! d( FIG(2) )/d(eta)
      SFD(2,3) = - 0.125D0 * (1.0D0 - xez(2)) * (1.0D0 + xez(1))     ! d( FIG(2) )/d(zeta)

      SFD(3,1) =   0.125D0 * (1.0D0 - xez(3)) * (1.0D0 + xez(2))     ! d( FIG(3) )/d(xi)
      SFD(3,2) =   0.125D0 * (1.0D0 - xez(3)) * (1.0D0 + xez(1))     ! d( FIG(3) )/d(eta)
      SFD(3,3) = - 0.125D0 * (1.0D0 + xez(2)) * (1.0D0 + xez(1))     ! d( FIG(3) )/d(zeta)

      SFD(4,1) = - 0.125D0 * (1.0D0 - xez(3)) * (1.0D0 + xez(2))     ! d( FIG(4) )/d(xi)
      SFD(4,2) =   0.125D0 * (1.0D0 - xez(3)) * (1.0D0 - xez(1))     ! d( FIG(4) )/d(eta)
      SFD(4,3) = - 0.125D0 * (1.0D0 + xez(2)) * (1.0D0 - xez(1))     ! d( FIG(4) )/d(zeta)

      SFD(5,1) = - 0.125D0 * (1.0D0 + xez(3)) * (1.0D0 - xez(2))     ! d( FIG(5) )/d(xi)
      SFD(5,2) = - 0.125D0 * (1.0D0 + xez(3)) * (1.0D0 - xez(1))     ! d( FIG(5) )/d(eta)
      SFD(5,3) =   0.125D0 * (1.0D0 - xez(2)) * (1.0D0 - xez(1))     ! d( FIG(5) )/d(zeta)

      SFD(6,1) =   0.125D0 * (1.0D0 + xez(3)) * (1.0D0 - xez(2))     ! d( FIG(6) )/d(xi)
      SFD(6,2) = - 0.125D0 * (1.0D0 + xez(3)) * (1.0D0 + xez(1))     ! d( FIG(6) )/d(eta)
      SFD(6,3) =   0.125D0 * (1.0D0 - xez(2)) * (1.0D0 + xez(1))     ! d( FIG(6) )/d(zeta)

      SFD(7,1) =   0.125D0 * (1.0D0 + xez(3)) * (1.0D0 + xez(2))     ! d( FIG(7) )/d(xi)
      SFD(7,2) =   0.125D0 * (1.0D0 + xez(3)) * (1.0D0 + xez(1))     ! d( FIG(7) )/d(eta)
      SFD(7,3) =   0.125D0 * (1.0D0 + xez(2)) * (1.0D0 + xez(1))     ! d( FIG(7) )/d(zeta)

      SFD(8,1) = - 0.125D0 * (1.0D0 + xez(3)) * (1.0D0 + xez(2))     ! d( FIG(8) )/d(xi)
      SFD(8,2) =   0.125D0 * (1.0D0 + xez(3)) * (1.0D0 - xez(1))     ! d( FIG(8) )/d(eta)
      SFD(8,3) =   0.125D0 * (1.0D0 + xez(2)) * (1.0D0 - xez(1))     ! d( FIG(8) )/d(zeta)


      END

! ------------------------------------------------------------------------------------    
      SUBROUTINE hex_kks2lks_org(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                        x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xi,eta,zeta,ierr)

!
!     $: na podlagi x,y,z najde xi,eta,zeta znotraj 8ih vogalov
!        glej NR stran 372 poglavje 9-6 Newton-Rapson method
!
!        ierr=0, tocka (x,y,z) je v heksaedru, ki ga doloca 8 vozlisc
!        ierr=1, tocka (x,y,z) NI v heksaedru, ki ga doloca 8 vozlisc
!
!          7 ---- 8
!         /|     /|
!        / |    / |
!       5 ---- 6  |
!       |  |   |  |
!       |  3 --|- 4
!       | /    |/
!       1 ---- 2
!      
! ------------------------------------------------------------------------------------    
      REAL(8) xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z
      REAL(8) xx(3),eps
      INTEGER maxit,nit,ierr

!     zacetni priblizek
      xx=0.0D0
!     toleranca 
      eps=1.0D-10 !-6 !-13
!     najvecje stevilo iteracij
      maxit=1000 !10000
            
      CALL mnewt(maxit,xx,3,eps,eps,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                              x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,nit)

!     zapomnimo rezultate, preverimo ce skonvergiralo !!
      IF (nit.GE.maxit) THEN
!        WRITE (*,*) "Dosezeno maximalno stevilo iteracij, kks2lks",nit
!        print *,xi,eta,zeta
         ierr=1
         RETURN
      END IF

      xi=xx(1)
      eta=xx(2)
      zeta=xx(3)       

!     preverimo, ce je noter ali zunaj
      IF ( xi.GE.-1.0D0.AND.  xi.LE.1.0D0.AND. &
     &    eta.GE.-1.0D0.AND. eta.LE.1.0D0.AND. &
        zeta.GE.-1.0D0.AND.zeta.LE.1.0D0) THEN
        ierr=0
      ELSE 
        ierr=1 
!        print *,xi,eta,zeta,nit
!        print *,"kk2lks",x,y,z
!        print *,x1,y1,z1
!        print *,x2,y2,z2
!        print *,x3,y3,z3
!        print *,x4,y4,z4
!        print *,x5,y5,z5
!        print *,x6,y6,z6
!        print *,x7,y7,z7
!        print *,x8,y8,z8
!        print *," "
      END IF

      END

! ------------------------------------------------------------------------------------    
      SUBROUTINE hex_kks2lks(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                         x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xi,eta,zeta,ierr)

!
!     $: na podlagi x,y,z najde xi,eta,zeta znotraj 8ih vogalov
!        glej NR stran 372 poglavje 9-6 Newton-Rapson method
!
!        ierr=0, tocka (x,y,z) je v heksaedru, ki ga doloca 8 vozlisc
!        ierr=1, tocka (x,y,z) NI v heksaedru, ki ga doloca 8 vozlisc
!
!          7 ---- 8
!         /|     /|
!        / |    / |
!       5 ---- 6  |
!       |  |   |  |
!       |  3 --|- 4
!       | /    |/
!       1 ---- 2
!      
! ------------------------------------------------------------------------------------    
      REAL(8) xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z
      REAL(8) xx(3),eps

      REAL(8) minVal, maxVal

      INTEGER maxit,nit,ierr

!     zacetni priblizek
      xx=0.0D0
!     toleranca 
      eps=1.0D-10 !-6 !-13
!     najvecje stevilo iteracij
      maxit=1000 !10000
            
      CALL mnewt(maxit,xx,3,eps,eps,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                              x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,nit)

!     zapomnimo rezultate, preverimo ce skonvergiralo !!
      IF (nit.GE.maxit) THEN
!        WRITE (*,*) "Dosezeno maximalno stevilo iteracij, kks2lks",nit
!        print *,xi,eta,zeta
         ierr=1
         RETURN
      END IF

      xi=xx(1)
      eta=xx(2)
      zeta=xx(3)   
      
      minVal = -1.0D0 - 5e-01
      maxVal = 1.0D0 + 5e-01

!     preverimo, ce je noter ali zunaj
      IF ( xi.GE.minVal.AND.  xi.LE.maxVal.AND. &
      eta.GE.minVal.AND. eta.LE.maxVal.AND. &
      zeta.GE.minVal.AND.zeta.LE.maxVal) THEN
        ierr=0
      ELSE 
        ierr=1 
!        print *,"error in hex_kks2lks"       
!        print *,xi,eta,zeta,nit
!        print *,"kk2lks",x,y,z
!        print *,x1,y1,z1
!        print *,x2,y2,z2
!        print *,x3,y3,z3
!        print *,x4,y4,z4
!        print *,x5,y5,z5
!        print *,x6,y6,z6
!        print *,x7,y7,z7
!        print *,x8,y8,z8
!        print *," "
      END IF

      END      
      
! ------------------------------------------------------------------------------------    
      SUBROUTINE kkk2lks_userfun(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                                x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xi,eta,zeta,F,J)

!
!     $: na podlagi x,y,z in priblizka xi,eta,zeta izracuna f, ki ga
!       minimiziramo in odvode jacobijeve 
!      
! ------------------------------------------------------------------------------------    
      REAL(8) xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z

      REAL(8) FIG(8),eta1m,eta1p,eta2m,eta2p,eta3m,eta3p
      
      REAL(8) F(3),J(3,3)

!      xi=1.0
!      eta=1.0
!      zeta=1.0

      ETA1M=1.0D0-xi
      ETA1P=1.0D0+xi   
      ETA2M=1.0D0-eta
      ETA2P=1.0D0+eta   
      ETA3M=1.0D0-zeta   
      ETA3P=1.0D0+zeta  

      FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
      FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
      FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
      FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
      FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
      FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
      FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
      FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

!      print *,FIg
!      stop

!     funkcije, ki jih minimiziramo
      F(1)=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4+FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8-X
      F(2)=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4+FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8-Y
      F(3)=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4+FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8-Z

!     Jacobijeva matrika J_ij=\p F_i / \p x_j      
      J(1,1)=-ETA2M*ETA3M*X1+ETA2M*ETA3M*X2-ETA2P*ETA3M*X3+ETA2P*ETA3M*X4 &
            -ETA2M*ETA3P*X5+ETA2M*ETA3P*X6-ETA2P*ETA3P*X7+ETA2P*ETA3P*X8

      J(2,1)=-ETA2M*ETA3M*Y1+ETA2M*ETA3M*Y2-ETA2P*ETA3M*Y3+ETA2P*ETA3M*Y4 &
            -ETA2M*ETA3P*Y5+ETA2M*ETA3P*Y6-ETA2P*ETA3P*Y7+ETA2P*ETA3P*Y8

      J(3,1)=-ETA2M*ETA3M*Z1+ETA2M*ETA3M*Z2-ETA2P*ETA3M*Z3+ETA2P*ETA3M*Z4 &
            -ETA2M*ETA3P*Z5+ETA2M*ETA3P*Z6-ETA2P*ETA3P*Z7+ETA2P*ETA3P*Z8
!     
      J(1,2)=-ETA1M*ETA3M*X1-ETA1P*ETA3M*X2+ETA1M*ETA3M*X3+ETA1P*ETA3M*X4 &
            -ETA1M*ETA3P*X5-ETA1P*ETA3P*X6+ETA1M*ETA3P*X7+ETA1P*ETA3P*X8

      J(2,2)=-ETA1M*ETA3M*Y1-ETA1P*ETA3M*Y2+ETA1M*ETA3M*Y3+ETA1P*ETA3M*Y4 &
            -ETA1M*ETA3P*Y5-ETA1P*ETA3P*Y6+ETA1M*ETA3P*Y7+ETA1P*ETA3P*Y8

      J(3,2)=-ETA1M*ETA3M*Z1-ETA1P*ETA3M*Z2+ETA1M*ETA3M*Z3+ETA1P*ETA3M*Z4 &
            -ETA1M*ETA3P*Z5-ETA1P*ETA3P*Z6+ETA1M*ETA3P*Z7+ETA1P*ETA3P*Z8
!
      J(1,3)=-ETA1M*ETA2M*X1-ETA1P*ETA2M*X2-ETA1M*ETA2P*X3-ETA1P*ETA2P*X4 &
            +ETA1M*ETA2M*X5+ETA1P*ETA2M*X6+ETA1M*ETA2P*X7+ETA1P*ETA2P*X8

      J(2,3)=-ETA1M*ETA2M*Y1-ETA1P*ETA2M*Y2-ETA1M*ETA2P*Y3-ETA1P*ETA2P*Y4 &
            +ETA1M*ETA2M*Y5+ETA1P*ETA2M*Y6+ETA1M*ETA2P*Y7+ETA1P*ETA2P*Y8

      J(3,3)=-ETA1M*ETA2M*Z1-ETA1P*ETA2M*Z2-ETA1M*ETA2P*Z3-ETA1P*ETA2P*Z4 &
            +ETA1M*ETA2M*Z5+ETA1P*ETA2M*Z6+ETA1M*ETA2P*Z7+ETA1P*ETA2P*Z8      
      
      END

! ------------------------------------------------------------------------------------
      SUBROUTINE hex_lks2kks(xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                                    x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z)
!
!     transformacija iz lokalnega koor. sistema (xi,eta,zeta) v kartezijevega (x,y,z)
!     za heksaeder dolocen z 8 vogali
!    
! ------------------------------------------------------------------------------------          
      REAL(8) xi,eta,zeta,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z
      REAL(8) FIG(8),eta1m,eta1p,eta2m,eta2p,eta3m,eta3p

      ETA1M=1.0D0-xi
      ETA1P=1.0D0+xi   
      ETA2M=1.0D0-eta
      ETA2P=1.0D0+eta   
      ETA3M=1.0D0-zeta   
      ETA3P=1.0D0+zeta  

      FIG(1)=0.125D0*ETA1M*ETA2M*ETA3M
      FIG(2)=0.125D0*ETA1P*ETA2M*ETA3M
      FIG(3)=0.125D0*ETA1M*ETA2P*ETA3M
      FIG(4)=0.125D0*ETA1P*ETA2P*ETA3M
      FIG(5)=0.125D0*ETA1M*ETA2M*ETA3P
      FIG(6)=0.125D0*ETA1P*ETA2M*ETA3P
      FIG(7)=0.125D0*ETA1M*ETA2P*ETA3P
      FIG(8)=0.125D0*ETA1P*ETA2P*ETA3P

      X=FIG(1)*X1+FIG(2)*X2+FIG(3)*X3+FIG(4)*X4+FIG(5)*X5+FIG(6)*X6+FIG(7)*X7+FIG(8)*X8
      Y=FIG(1)*Y1+FIG(2)*Y2+FIG(3)*Y3+FIG(4)*Y4+FIG(5)*Y5+FIG(6)*Y6+FIG(7)*Y7+FIG(8)*Y8
      Z=FIG(1)*Z1+FIG(2)*Z2+FIG(3)*Z3+FIG(4)*Z4+FIG(5)*Z5+FIG(6)*Z6+FIG(7)*Z7+FIG(8)*Z8
      
      END
     


      
! ------------------------------------------------------------------------------------        
      SUBROUTINE mnewt(ntrial,xx,n,tolx,tolf,x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                              x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,nit)
      INTEGER n,ntrial
      REAL(8) tolf,tolx,xx(n)
      REAL(8) x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4
      REAL(8) x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,x,y,z      
!U    USES lubksb,ludcmp,usrfun
      INTEGER i,k,indx(n),nit
      REAL(8) d,errf,errx,fjac(n,n),fvec(n),p(n)
      do 14  k=1,ntrial
        nit=k
        CALL kkk2lks_userfun(x,y,z,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4, &
                     x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,xx(1),xx(2),xx(3),Fvec,Fjac)
        errf=0.
        do 11 i=1,n
          errf=errf+abs(fvec(i))
11      continue
        if(errf.le.tolf)return
        do 12 i=1,n
          p(i)=-fvec(i)
12      continue
        call ludcmp(fjac,n,n,indx,d)
        call lubksb(fjac,n,n,indx,p)
        errx=0.
        do 13 i=1,n
          errx=errx+abs(p(i))
          xx(i)=xx(i)+p(i)
13      continue
        if(errx.le.tolx)return
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.  
! ------------------------------------------------------------------------------------        
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER n,np,indx(n)
      REAL(8) a(np,np),b(n)
      INTEGER i,ii,j,ll
      REAL(8) sum
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.
! ------------------------------------------------------------------------------------      
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER n,np,indx(n),NMAX
      REAL(8) d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)
      INTEGER i,imax,j,k
      REAL(8) aamax,dum,sum,vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) then
          print *,'singular matrix in ludcmp'
          stop
        end if
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software *5ji6.)+.
