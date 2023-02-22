SUBROUTINE probe(v,itc)
    USE COMUM
    IMPLICIT NONE
    INTEGER :: i,j,itc

    REAL(8), DIMENSION(1:imax  ,1:jmax) :: v
    REAL(8) :: st_pos

    DO j=1,jmax
         IF (v(1,j) .LT. 0.d0 .AND. v(1,j+1) .GT. 0.d0) THEN
             st_pos = 0.5d0 * (y(j) + y(j+1))
         ENDIF
    ENDDO

    open(unit=550,file='data/probe.dat',status='unknown',position='append')

        write(550,*) itc, st_pos-0.5d0

    close(550)



RETURN
END SUBROUTINE probe



