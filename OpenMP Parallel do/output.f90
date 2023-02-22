SUBROUTINE output(um,vm,u,v,p,Z,T,H,k)
    USE COMUM
    use omp_lib
    IMPLICIT NONE
    INTEGER :: i,j,k

    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: u,v,P,Z,T,H,Yi
    character*16 filename
    character(len=100) :: filename2 
!!$omp parallel do private(i, j)
    do i=1,imax
        do j=1,jmax
     
        if(Z(i,j).gt.(1.d0)) then
            Yi(i,j) = (Z(i,j) - 1.d0) / S 
        else
            Yi(i,j) = - Z(i,j) + 1.d0 
        endif

        enddo
     enddo
!!$omp end parallel do


    open (550,file='data/output_variables.dat')
        
         do i = 1,imax
             do j = 1,jmax
            
             write (550,*) X(i) , Y(j) , u(i,j), v(i,j) , p(i,j) ,&
                             T(i,j)*Too , Z(i,j) , H(i,j), Yi(i,j)

             enddo
         enddo    

    close(550)


    !############ RESTART ###############

    write(filename,*) k
    filename2 = adjustl(filename)    
    open (550,file='data/results/U/restartU'//trim(filename2)//'.dat',status='unknown')


        do i=1,imax+1
            do j = 1,jmax
        
            write (550,*) um(i,j)
        
            enddo
        enddo

    close(550)

    open (550,file='data/results/V/restartV'//trim(filename2)//'.dat',status='unknown')


         do i=1,imax
            do j = 1,jmax+1
        
            write (550,*) vm(i,j)
        
            enddo
         enddo

    close(550)

    open (550,file='data/results/PTZH/restartPTZH'//trim(filename2)//'.dat',status='unknown')


        do i=1,imax
            do j = 1,jmax
        
            write (550,*) p(i,j),T(i,j),Z(i,j),H(i,j)
        
            enddo
        enddo

    close(550)

    !########## RESTART/RESTART.dat

    open (550,file='data/restart/restartU.dat',status='unknown')


        do i=1,imax+1
            do j = 1,jmax
        
            write (550,*) um(i,j)
        
            enddo
        enddo

    close(550)

    open (550,file='data/restart/restartV.dat',status='unknown')


         do i=1,imax
            do j = 1,jmax+1
        
            write (550,*) vm(i,j)
        
            enddo
         enddo

    close(550)

    open (550,file='data/restart/restartPTZH.dat',status='unknown')


        do i=1,imax
            do j = 1,jmax
        
            write (550,*) p(i,j),T(i,j),Z(i,j),H(i,j)
        
            enddo
        enddo

    close(550)

RETURN
END SUBROUTINE output



