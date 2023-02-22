PROGRAM pos_process
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: u,v,P,Z,T,H,Yi
    character*16 filename
    character(len=100) :: filename2

    CALL properties

!-------- reference values --------------   
     OPEN(unit=550,file='input/reference.dat',status='old')
        READ(550,nml=ref)
        WRITE(*,nml=ref)
    CLOSE(unit=550)


!############ READ THE INPUT DATA FROM OUTPUT.F90 ##############################
    OPEN (550,file='data/output_variables.dat',status='old',access='sequential')

         do i = 1,imax
             do j = 1,jmax
            
             read (550,*) X(i) , Y(j) , u(i,j), v(i,j) , p(i,j) ,&
                             T(i,j) , Z(i,j) , H(i,j), Yi(i,j)

             enddo
         enddo    

     CLOSE(550)


    OPEN (550,file='data/mesh.dat')

        WRITE (550,*) imax, jmax, 10

    CLOSE(550)

!########## OUTPUT TO PLOT FIELDS IN PYTHON ################################
    open(unit=550, file='data/x.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)  x(i)
          end do
    CLOSE(550)
    
    
    open(unit=550, file='data/y.dat', ACTION="write", STATUS="replace")
          do j=1,jmax
          write(550, *)  y(j)
          end do
    CLOSE(550)
    
    
    open(unit=550, file='data/u.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(u(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)

    open(unit=550, file='data/v.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(v(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)
    
    open(unit=550, file='data/P.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(P(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)
    
    open(unit=550, file='data/T.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(T(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)
    
    open(unit=550, file='data/Z.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(Z(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)



!########## OUTPUT TO USE IN PARAVIEW APP ################################
    OPEN(unit=1,file='data/paraview_output.vtk',status='unknown')

    WRITE(1,'(a)')'# vtk DataFile Version 3.0' 
    WRITE(1,'(a)')'Droplet Combustion'   
    WRITE(1,'(a)')'ASCII'   
    WRITE(1,'(a)')'DATASET STRUCTURED_GRID'  
    WRITE(1,'(A10,A1,I3,A1,I3,2A1)')'DIMENSIONS',' ',jmax,' ',imax,' ','1'
!WRITE(1,fmt='(2i4.2)') nj,ni

!WRITE(1,'(A10,I2,I2,A1)')'DIMENSIONS', nj,ni,'1'
!WRITE(1,'(A,I8.0)')'DIMENSIONS', nj

    WRITE(1,'(A6,A1,I6,A1,A5)')'POINTS',' ', jmax*imax,' ','float'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(3F14.4)') x(i), y(j), 0.d0
        ENDDO
    ENDDO 
    WRITE(1,'(A10,A1,I6)')'POINT_DATA',' ', imax*jmax
    WRITE(1,'(a)')'SCALARS VEL_MAGNITUDE float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') dsqrt((u(i,j))**2+(v(i,j))**2)
        ENDDO
    ENDDO

    WRITE(1,'(a)')'SCALARS T float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') 	T(i,j)
        ENDDO
    ENDDO

    WRITE(1,'(a)')'SCALARS P float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') 	P(i,j)
        ENDDO
    ENDDO

    WRITE(1,'(a)')'SCALARS Z float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') 	Z(i,j)
        ENDDO
    ENDDO

    WRITE(1,'(a)')'SCALARS H float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') 	H(i,j)
        ENDDO
    ENDDO

    WRITE(1,'(a)')'SCALARS Yi float'
    WRITE(1,'(a)')'LOOKUP_TABLE default'
    DO i=1,imax
        DO j=1,jmax
            WRITE(1,'(F14.4)') 	Yi(i,j)
        ENDDO
    ENDDO

    WRITE(1,'(a)')''
    WRITE(1,'(a)')''
    WRITE(1,'(a)')'VECTORS Vectors float'

    DO i=1,imax
        DO j=1,jmax		
            WRITE(1,'(3F14.4)') 	u(i,j), v(i,j), 0.d0
        ENDDO
    ENDDO

    CLOSE(1)


RETURN
END PROGRAM pos_process
