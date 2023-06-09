subroutine transient(u,v,p,T,Z,tr)
    use comum
    implicit none
    integer :: i, j, tr
    REAL(8), DIMENSION(1:imax,1:jmax) :: u, v, P, Z, T
    
    character(len=100) :: filename2

    write(filename2,*) (tr/n_tr)+100
    filename2 = adjustl(filename2)   

    
    open(unit=550, file='transient/data/u'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(u(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)

    open(unit=550, file='transient/data/v'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(v(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)
    
    open(unit=550, file='transient/data/P'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(P(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)
    
    open(unit=550, file='transient/data/T'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(T(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)
    
    open(unit=550, file='transient/data/Z'//trim(filename2)//'.dat', ACTION="write", STATUS="replace")
          do i=1,imax
          write(550, *)( real(Z(i,j)) ,j=1,jmax)
          end do
    CLOSE(550)

    open (550,file='transient/data/time'//trim(filename2)//'.dat')

                write (550,*) time

    close(550)
    
    open (550,file='transient/data/grid.dat')


         do j=1,jmax
                do i = 1,imax
            
                write (550,*) X(i) , Y(j)

                enddo
         enddo    

    close(550)

!    i = 1
!    do j=int(y_down/dy),jmax
!                
!        if (Zc(i,j) .gt. 1.d0) then
!
!        yf = ((1.d0-Zc(i,j-1))*( Y(j)- Y(j-1) )) / ( Zc(i,j) - Zc(i,j-1) ) + Y(j-1)
!        
!        go to 200        
!
!        endif    
!        
!    200 enddo    
!
!    open(unit=550,file='data/flametip.dat',status='unknown',position='append')
!        
!        if (yf .lt. y_up) then
!            write(550,*) itc,yf
!        else
!            write(550,*) itc, y_up
!        endif
!              
!        
!    close(550)
!
return
end subroutine transient


