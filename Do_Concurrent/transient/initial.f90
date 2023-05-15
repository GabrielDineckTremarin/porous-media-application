SUBROUTINE init
    USE comum
    IMPLICIT NONE
    
!-------- iterations values -------------    
    OPEN(unit=550,file='input/iterations.dat',status='old')
        READ(550,nml=iterations)
        WRITE(*,nml=iterations)
    CLOSE(unit=550)
    
    
!-------- reference values --------------   
     OPEN(unit=550,file='input/reference.dat',status='old')
        READ(550,nml=ref)
        WRITE(*,nml=ref)
    CLOSE(unit=550)
    
    
END SUBROUTINE init

SUBROUTINE IC(um,vm,p,T)
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j

    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: P, T
    REAL(8) :: eta_e

    vm = 0.1d0 !*v_i
    um = 0.d0
    p = 1.0d0
    T = Tinf
   ! Z = 0.d0


!    DO i=1,imax
!        P(i,1:jmax) = -x(i)*(1.d0) + Lhori + 1.d0  
!    ENDDO

!IC condition for stokes second problem
!    DO j=1,jmax
!        eta_e = dsqrt(6.0*Re/2.d0)*y(j)
!        um(1:imax,j) = dexp(-eta_e) * dcos(-eta_e)   
!    ENDDO

!    DO j=1,jmax
!            DO i = 1,imax
!
!            IF (Y(j) .GT. 0.d0) THEN
!
!            Z(i,j) = ( (S+1.d0) / dsqrt(Y(j)) ) &
!                       * dexp( - ( 1.d0 * v_i * Pe / 4.d0 ) * X(i)**2.d0 / (Y(j)) )
!
!            ENDIF
!
!            ENDDO
!    ENDDO

    call system ('rm data/flametip.dat')
    call system ('rm data/error.dat')

RETURN
END SUBROUTINE IC


SUBROUTINE restart(um,vm,p,T)
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T,Yi

    WRITE(*,*) 'RESTARTING PROGRAM'

    OPEN (550,file='data/restart/restartU.dat',status='old',access='sequential')

        do i=1,imax+1
            do j = 1,jmax
        
            read (550,*) um(i,j)
        
            enddo
        enddo

    CLOSE(550)

    OPEN (550,file='data/restart/restartV.dat',status='old',access='sequential')

         do i=1,imax
            do j = 1,jmax+1
        
            read (550,*) vm(i,j)
        
            enddo
         enddo

    CLOSE(550)

    OPEN (550,file='data/restart/restartPTZH.dat',status='old',access='sequential')

        do i=1,imax
            do j = 1,jmax
        
            read (550,*) p(i,j),T(i,j)
        
           enddo
        enddo

    CLOSE(550)


RETURN
END SUBROUTINE restart



SUBROUTINE restart_dom(um,vm,p,T,Z,H)
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER, parameter :: rimax = 41 !em x
    INTEGER, parameter :: rjmax = 321 !em y
    REAL(8), DIMENSION(1:rimax+1,1:rjmax  ) :: umr
    REAL(8), DIMENSION(1:rimax  ,1:rjmax+1) :: vmr
    REAL(8), DIMENSION(1:rimax,1:rjmax) :: Pres,Zr,Tr,Hr,H_res
    
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,Z,T,H


    WRITE(*,*) 'RESTARTING PROGRAM'

    OPEN (550,file='data/restart/restartU.dat',status='old',access='sequential')

        do i=1,rimax+1
            do j = 1,rjmax
        
            read (550,*) umr(i,j)
        
            enddo
        enddo

    CLOSE(550)

    OPEN (550,file='data/restart/restartV.dat',status='old',access='sequential')

         do i=1,rimax
            do j = 1,rjmax+1
        
            read (550,*) vmr(i,j)
        
            enddo
         enddo

    CLOSE(550)

    OPEN (550,file='data/restart/restartPTZH.dat',status='old',access='sequential')

        do i=1,rimax
            do j = 1,rjmax
        
            read (550,*) Pres(i,j),Tr(i,j),Zr(i,j),H_res(i,j)
        
            Hr(i,j) = H_res(i,j) + ( ((S + 1.d0) * Lf * Tinf / q + 1.d0) - H_res(i,j) )

            enddo
        enddo

    CLOSE(550)

    DO j=1,jmax
            DO i = 1,imax

                IF (j .LE. rjmax) THEN
       
                    p(i,j) = pres(i,j)
                    T(i,j) = Tr(i,j)
                    Z(i,j) = Zr(i,j)
                    H(i,j) = Hr(i,j)
       
                ELSE
       
                    p(i,j) = p(i,j-1)
                    T(i,j) = T(i,j-1)
                    Z(i,j) = Z(i,j-1)
                    H(i,j) = H(i,j-1)
       
                ENDIF
           ENDDO
    ENDDO


    DO j=1,jmax
            DO i = 1,imax+1

                IF (j .LE. rjmax) THEN
       
                    um(i,j) = umr(i,j)
       
                ELSE
       
                    um(i,j) = um(i,j-1)
       
                ENDIF
           ENDDO
    ENDDO


    DO j=1,jmax+1
            DO i = 1,imax

                IF (j .LE. rjmax) THEN
       
                    vm(i,j) = vmr(i,j)
       
                ELSE
       
                    vm(i,j) = vm(i,j-1)
       
                ENDIF
           ENDDO
    ENDDO


RETURN
END SUBROUTINE restart_dom
