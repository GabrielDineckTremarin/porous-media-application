
#line 1 "/home/gabriel/POROUS-MEDIA-APPLICATION-main/OpenMP_Teams/flametip.f90"
SUBROUTINE flametip(Z,itc)
    USE COMUM
    IMPLICIT NONE
      include 'flametip.f90.opari.inc'
#line 4 "/home/gabriel/POROUS-MEDIA-APPLICATION-main/OpenMP_Teams/flametip.f90"
    INTEGER :: i,j,itc

    REAL(8), DIMENSION(1:imax  ,1:jmax) :: Z
    REAL(8) :: yf


    i = 1
    do j=int(y_down/dx_c),jmax
                
        if (Z(i,j) .gt. 1.d0) then

        yf = ((1.d0-Z(i,j-1))*( Y(j)- Y(j-1) )) / ( Z(i,j) - Z(i,j-1) ) + Y(j-1)
        
        go to 200        

        endif    
        
    200 enddo    

    open(unit=550,file='data/flametip.dat',status='unknown',position='append')
        
        if (yf .lt. y_up) then
            write(550,*) itc,yf
        else
            write(550,*) itc, y_up
        endif
              
        
    close(550)


RETURN
END SUBROUTINE flametip



