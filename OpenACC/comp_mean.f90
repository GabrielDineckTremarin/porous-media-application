SUBROUTINE comp_mean(u,v,um,vm)
    USE comum
    ! $acc routine
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: u, v

    ! $acc data present(um, vm, u, v)

!##############CALCULA OS PONTOS MEDIOS DAS VELOCIDADES ######################
!!$omp parallel sections private (i, j)
!!$omp section
! $acc parallel loop collapse(2)
    DO i=1,imax
        DO j=1,jmax
            u(i,j) = (um(i+1,j)+um(i,j))*0.5d0
        ENDDO
    ENDDO
!!$omp section

! $acc parallel loop collapse(2)
    DO i=1,imax
        DO j=1,jmax
            v(i,j) = (vm(i,j+1)+vm(i,j))*0.5d0
        eNDDO
    ENDDO
!!$omp end parallel sections

! $acc end data
RETURN
END SUBROUTINE comp_mean
