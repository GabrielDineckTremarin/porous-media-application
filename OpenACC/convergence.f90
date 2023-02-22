!!! check convergence
SUBROUTINE convergence(itc,c2,error,residual_p,residual_u,residual_v)
    use comum
    use openacc
    !acc routine
    IMPLICIT NONE
    
    integer :: i, j, itc
    real(8) :: c2, error

    real(8) :: residual_p,residual_u,residual_v

!$acc data present(itc, error, residual_u, residual_v, residual_p)
    itc = itc+1

    error = MAX(residual_u,residual_v,residual_p)
!$acc end data

    !IF (MOD(itc,100).EQ.0) THEN

    !    WRITE(*,*) itc,' ',residual_u,residual_v,residual_p
       
    !ENDIF

    !    OPEN(unit=549,file='data/error.dat',status='unknown',position='appEND')
    !         WRITE(549,*) itc,' ',residual_u,residual_v,residual_p
    !    close(549)
RETURN
END SUBROUTINE convergence

SUBROUTINE iterat_convergence(itc, error)
    use comum
    use openacc
    !acc routine
    IMPLICIT NONE

    integer :: itc
    real(8) ::  error

!$acc data present(itc, error)
    itc = 0
    error = 100.d0
!$acc end data

    !IF (MOD(itc,100).EQ.0) THEN

    !    WRITE(*,*) itc,' ',residual_u,residual_v,residual_p

    !ENDIF

    !    OPEN(unit=549,file='data/error.dat',status='unknown',position='appEND')
    !         WRITE(549,*) itc,' ',residual_u,residual_v,residual_p
    !    close(549)
RETURN
END SUBROUTINE iterat_convergence
