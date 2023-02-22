SUBROUTINE comp_mean(u,v,um,vm)
    USE comum
    !use omp_lib
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: u, v


!##############CALCULA OS PONTOS MEDIOS DAS VELOCIDADES ######################
!!$omp parallel sections private (i, j)
!!$omp section
    DO i=1,imax
        DO j=1,jmax
            u(i,j) = (um(i+1,j)+um(i,j))*0.5d0
        ENDDO
    ENDDO
!!$omp section
    DO i=1,imax
        DO j=1,jmax
            v(i,j) = (vm(i,j+1)+vm(i,j))*0.5d0
        eNDDO
    ENDDO
!!$omp end parallel sections

RETURN
END SUBROUTINE comp_mean

