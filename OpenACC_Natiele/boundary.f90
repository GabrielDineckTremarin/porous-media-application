subroutine bcUV(um,vm)
    use comum
    use omp_lib
    implicit none
    integer :: i, j
    real(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    real(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
       
!------ contorno inferior -----------
    !$omp parallel sections private (i)
    !$omp section  
    DO i=1,imax
        vm(i,2) = v_i
        vm(i,1) = v_i
    ENDDO
    !$omp section  
    DO i=2,imax
        um(i,1) = 0.d0    !!!DARCY
    ENDDO
   !$omp end parallel sections

!------ contorno superior ----------
!$omp parallel sections private (i)
!$omp section 
    DO i=2,imax
        vm(i,jmax  ) = 2.d0*vm(i,jmax-2) - vm(i,jmax-1)
        vm(i,jmax+1) = 2.d0*vm(i,jmax-1) - vm(i,jmax)
    ENDDO
!$omp section 
    DO i=2,imax
        um(i,jmax) = um(i,jmax-1) 
    ENDDO
!$omp end parallel sections    
   

 !------ contorno esquerdo --------
!$omp parallel sections private (j)
!$omp section 
    DO j=1,jmax
        !um(2,j) = 0.d0 !symmetry
        !um(1,j) = 0.d0 !symmetry
        um(2,j) = 2.d0*um(4  ,j) - um(3,j)  !Darcy
        um(1,j) = 2.d0*um(3  ,j) - um(2,j)	!Darcy
      !    um(2,j) = 0.d0
      !  um(1,j) = 0.d0
    ENDDO
!$omp section 
    DO j=1,jmax+1
        vm(1,j) = vm(2,j) !darcy
    ENDDO
 !$omp end parallel sections

!------ contorno direito --------
!$omp parallel sections private (j)
!$omp section 
    DO j=1,jmax
        um(imax  ,j) = 2.d0*um(imax-2  ,j) - um(imax-1,j)
        um(imax+1,j) = 2.d0*um(imax-1  ,j) - um(imax,j)
    ENDDO
!$omp section 
    DO j=2,jmax
        vm(imax,j) =  vm(imax-1,j) !Darcy
    ENDDO            
 !$omp end parallel sections

return
end subroutine bcUV
!------------------------------------------------------------------------------

subroutine bcUV1(ui,vm_tau)
    use comum
    use openacc
    !acc routine
    implicit none
    integer :: i, j
    real(8), DIMENSION(1:imax+1,1:jmax  ) :: ui
    real(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
! $acc data present(ui, vm_tau)
!------ contorno inferior -----------
    !$acc parallel loop private (i)
    DO i=1,imax
        vm_tau(i,2) = v_i
        vm_tau(i,1) = v_i
    ENDDO
    !$acc end parallel loop
    !$acc parallel loop private(i)
    DO i=2,imax
        ui(i,1) = 0.d0    !!!DARCY
    ENDDO
   !$acc end parallel loop

!------ contorno superior ----------
!$acc parallel loop private (i)
    DO i=2,imax
        vm_tau(i,jmax  ) = 2.d0*vm_tau(i,jmax-2) - vm_tau(i,jmax-1)
        vm_tau(i,jmax+1) = 2.d0*vm_tau(i,jmax-1) - vm_tau(i,jmax)
    ENDDO
!$acc end parallel loop
!$acc parallel loop 
    DO i=2,imax
        ui(i,jmax) = ui(i,jmax-1)
    ENDDO
!$acc end parallel   

!------ contorno esquerdo --------
!$acc parallel loop 
    DO j=1,jmax
        ui(2,j) = 2.d0*ui(4  ,j) - ui(3,j)  !Darcy
        ui(1,j) = 2.d0*ui(3  ,j) - ui(2,j)      !Darcy
    ENDDO
!$acc end parallel loop
!$acc parallel loop 
    DO j=1,jmax+1
        vm_tau(1,j) = vm_tau(2,j) !darcy
    ENDDO
!$acc end parallel loop

!------ contorno direito --------
!$acc parallel loop 
    DO j=1,jmax
        ui(imax  ,j) = 2.d0*ui(imax-2  ,j) - ui(imax-1,j)
        ui(imax+1,j) = 2.d0*ui(imax-1  ,j) - ui(imax,j)
    ENDDO
!$acc end parallel loop
!$acc parallel loop  
    DO j=2,jmax
        vm_tau(imax,j) =  vm_tau(imax-1,j) !Darcy
    ENDDO
!$acc end parallel loop 
! $acc end data
return
end subroutine bcUV1

!-----------------------------------------------------
subroutine bcUV2(um_n_tau,vm_tau)
    use comum
    use openacc
    !acc routine
    implicit none
    integer :: i, j
    real(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n_tau
    real(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
! $acc data present(um_n_tau, vm_tau)
!------ contorno inferior -----------
    !$acc parallel loop private (i)
    DO i=1,imax
        vm_tau(i,2) = v_i
        vm_tau(i,1) = v_i
    ENDDO
    !$acc end parallel loop
    !$acc parallel loop private(i)
    DO i=2,imax
        um_n_tau(i,1) = 0.d0    !!!DARCY
    ENDDO
   !$acc end parallel loop

!------ contorno superior ----------
!$acc parallel loop private (i)
    DO i=2,imax
        vm_tau(i,jmax  ) = 2.d0*vm_tau(i,jmax-2) - vm_tau(i,jmax-1)
        vm_tau(i,jmax+1) = 2.d0*vm_tau(i,jmax-1) - vm_tau(i,jmax)
    ENDDO
!$acc end parallel loop
!$acc parallel loop 
    DO i=2,imax
        um_n_tau(i,jmax) = um_n_tau(i,jmax-1)
    ENDDO
!$acc end parallel   

!------ contorno esquerdo --------
!$acc parallel loop 
    DO j=1,jmax
        um_n_tau(2,j) = 2.d0*um_n_tau(4  ,j) - um_n_tau(3,j)  !Darcy
        um_n_tau(1,j) = 2.d0*um_n_tau(3  ,j) - um_n_tau(2,j)      !Darcy
    ENDDO
!$acc end parallel loop
!$acc parallel loop 
    DO j=1,jmax+1
        vm_tau(1,j) = vm_tau(2,j) !darcy
    ENDDO
!$acc end parallel loop

!------ contorno direito --------
!$acc parallel loop 
    DO j=1,jmax
        um_n_tau(imax  ,j) = 2.d0*um_n_tau(imax-2  ,j) - um_n_tau(imax-1,j)
        um_n_tau(imax+1,j) = 2.d0*um_n_tau(imax-1  ,j) - um_n_tau(imax,j)
    ENDDO
!$acc end parallel loop
!$acc parallel loop  
    DO j=2,jmax
        vm_tau(imax,j) =  vm_tau(imax-1,j) !Darcy
    ENDDO
!$acc end parallel loop 
! $acc end data
return
end subroutine bcUV2

!-----------------------------------------------------
!------------------------------------------------------
subroutine bcUV3(um_tau,vi)
    use comum
    use openacc
    !$acc routine
    implicit none
    integer :: i, j
    real(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    real(8), DIMENSION(1:imax  ,1:jmax+1) :: vi

!------ contorno inferior -----------
    ! $acc parallel loop 
    DO i=1,imax
        vi(i,2) = v_i
        vi(i,1) = v_i
    ENDDO
    ! $acc end parallel loop
    ! $acc parallel loop 
    DO i=2,imax
        um_tau(i,1) = 0.d0    !!!DARCY
    ENDDO
   ! $ acc end parallel loop

!------ contorno superior ----------
! $acc parallel loop
    DO i=2,imax
        vi(i,jmax  ) = 2.d0*vi(i,jmax-2) - vi(i,jmax-1)
        vi(i,jmax+1) = 2.d0*vi(i,jmax-1) - vi(i,jmax)
    ENDDO
! $acc end parallel loop
! $acc parallel loop
    DO i=2,imax
        um_tau(i,jmax) = um_tau(i,jmax-1)
    ENDDO
! $acc end parallel loop
  
!------ contorno esquerdo --------
! $acc parallel loop 
    DO j=1,jmax
        um_tau(2,j) = 2.d0*um_tau(4  ,j) - um_tau(3,j)  !Darcy
        um_tau(1,j) = 2.d0*um_tau(3  ,j) - um_tau(2,j)            
    ENDDO
! $acc end parallel loop
! $acc parallel loop
    DO j=1,jmax+1
        vi(1,j) = vi(2,j) !darcy
    ENDDO
! $acc end parallel loop

!------ contorno direito --------
! $acc parallel loop 
    DO j=1,jmax
        um_tau(imax  ,j) = 2.d0*um_tau(imax-2  ,j) - um_tau(imax-1,j)
        um_tau(imax+1,j) = 2.d0*um_tau(imax-1  ,j) - um_tau(imax,j)
    ENDDO
! $acc end parallel loop

! $acc parallel loop
    DO j=2,jmax
        vi(imax,j) =  vi(imax-1,j) !Darcy
    ENDDO
! $acc end parallel loop

return
end subroutine bcUV3

!------------------------------------------------------
subroutine bcUV4(um_tau,vm_n_tau)
    use comum
    use omp_lib
    implicit none
    integer :: i, j
    real(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    real(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n_tau

!------ contorno inferior -----------
    !$omp parallel sections private (i)
    !$omp section  
    DO i=1,imax
        vm_n_tau(i,2) = v_i
        vm_n_tau(i,1) = v_i
    ENDDO
    !$omp section  
    DO i=2,imax
        um_tau(i,1) = 0.d0    !!!DARCY
    ENDDO
   !$omp end parallel sections

!------ contorno superior ----------
!$omp parallel sections private (i)
!$omp section 
    DO i=2,imax
        vm_n_tau(i,jmax  ) = 2.d0*vm_n_tau(i,jmax-2) - vm_n_tau(i,jmax-1)
        vm_n_tau(i,jmax+1) = 2.d0*vm_n_tau(i,jmax-1) - vm_n_tau(i,jmax)
    ENDDO
!$omp section 
    DO i=2,imax
        um_tau(i,jmax) = um_tau(i,jmax-1)
    ENDDO
!$omp end parallel sections  
!------ contorno esquerdo --------
!$omp parallel sections private (j)
!$omp section 
    DO j=1,jmax
        um_tau(2,j) = 2.d0*um_tau(4  ,j) - um_tau(3,j)  !Darcy
        um_tau(1,j) = 2.d0*um_tau(3  ,j) - um_tau(2,j)            
    ENDDO
!$omp section 
    DO j=1,jmax+1
        vm_n_tau(1,j) = vm_n_tau(2,j) !darcy
    ENDDO
 !$omp end parallel sections

!------ contorno direito --------
!$omp parallel sections private (j)
!$omp section 
    DO j=1,jmax
        um_tau(imax  ,j) = 2.d0*um_tau(imax-2  ,j) - um_tau(imax-1,j)
        um_tau(imax+1,j) = 2.d0*um_tau(imax-1  ,j) - um_tau(imax,j)
    ENDDO
!$omp section 
    DO j=2,jmax
        vm_n_tau(imax,j) =  vm_n_tau(imax-1,j) !Darcy
    ENDDO
 !$omp end parallel sections

return
end subroutine bcUV4

!------------------------------------------------------

SUBROUTINE bcP(pn)
    USE comum
    use openacc
    ! $acc routine
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: pn
!$acc data present (pn)

    !!! boundary condition
!--------contorno inferior e superior --------
   !$acc parallel loop private(i) 
   do i=1,imax
        pn(i,1  ) = pn(i,2)
        pn(i,jmax) = pn(i,jmax-1) + 1.d0*(pn(i,jmax-1)-pn(i,jmax-2))
    enddo
!$acc end parallel loop
!-------contorno esquerdo e direito -------
!$acc parallel loop private(j) 
    do j=1,jmax
        pn(1  ,j) =  pn(2,j)
        pn(imax,j) = pn(imax-1,j)
    enddo
!$acc end parallel loop
!$acc end data
RETURN
END SUBROUTINE bcP
!------------------------------------------------------------------------------------------
SUBROUTINE bcP1(pi)
    USE comum
    use openacc
    ! $acc routine
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: pi
!$acc data present (pi)
    !!! boundary condition
!--------contorno inferior e superior --------
!$acc parallel loop private(i) 
    do i=1,imax
        pi(i,1  ) = pi(i,2)
        pi(i,jmax) = pi(i,jmax-1) + 1.d0*(pi(i,jmax-1)-pi(i,jmax-2))
    enddo
!$acc end parallel loop
!-------contorno esquerdo e direito -------
!$acc parallel loop private(j) 
    do j=1,jmax
        pi(1  ,j) =  pi(2,j)
        pi(imax,j) = pi(imax-1,j)
    enddo
!$acc end parallel loop
!$acc end data
RETURN
END SUBROUTINE bcP1

!---------------------------------------------------------------------------------------------
SUBROUTINE bcZ(Zn)
    USE comum
    use openacc
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: Zn,gradZ_x,gradZ_y

    !!! boundary condition

!--------contorno inferior e superior --------
! $acc parallel loop
    do i=1,imax
        Zn(i,1  ) = Zn(i,2)
        Zn(i,jmax) = Zn(i,jmax-1) + 1.d0*(Zn(i,jmax-1)-Zn(i,jmax-2))
    enddo
! $acc end parallel loop
!-------contorno esquerdo e direito -------
! $acc parallel loop    
    do j=1,jmax
        Zn(1  ,j) = Zn(2,j)
        Zn(imax,j) = Zn(imax-1,j)
    enddo
! $acc end parallel loop
       i=1
! $acc parallel loop
            do j=2,jmax-1
            gradZ_x(i,j) = x(i) * ( Zn(i+1,j) - Zn(i,j  ) ) * dx(i)
            gradZ_y(i,j) = y(j) * ( Zn(i,j+1) - Zn(i,j  ) ) * dy(j)
            enddo
! $acc end parallel loop
       
! $acc parallel loop        
        do i=1,imax-1
             j=1
            gradZ_x(i,j) = x(i) * ( Zn(i+1,j) - Zn(i,j  ) ) * dx(i+1)
            gradZ_y(i,j) = y(i) * ( Zn(i,j+1) - Zn(i,j  ) ) * dy(j+1)
        enddo
! $acc end parallel loop
! $acc parallel loop collapse(2) private(i, j)
         do i=2,imax
            do j=2,jmax
            gradZ_x(i,j) = x(i) * ( Zn(i,j) - Zn(i-1,j  ) ) * dx(i)
            gradZ_y(i,j) = y(j) * ( Zn(i,j) - Zn(i  ,j-1) ) * dy(j)
            enddo
        enddo
! $acc END PARALLEL loop
! $acc end data        
RETURN
END SUBROUTINE bcZ

!------------------------------------------------------------------------------------
SUBROUTINE bcZ2(Zi)
    USE comum
    use openacc
! $acc routine
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: Zi,gradZ_x,gradZ_y
!$acc data present (Zi,x, y,dx,dy)
    !!! boundary condition
!--------contorno inferior e superior --------
!$acc parallel loop
    do i=1,imax
        Zi(i,1  ) = Zi(i,2)
        Zi(i,jmax) = Zi(i,jmax-1) + 1.d0*(Zi(i,jmax-1)-Zi(i,jmax-2))
    enddo
!$acc end parallel
!-------contorno esquerdo e direito -------
!$acc parallel loop
    do j=1,jmax
        Zi(1  ,j) = Zi(2,j)
        Zi(imax,j) = Zi(imax-1,j)
    enddo
!$acc end parallel loop
       i=1
!$acc parallel loop
            do j=2,jmax-1
            gradZ_x(i,j) = x(i) * ( Zi(i+1,j) - Zi(i,j  ) ) * dx(i)
            gradZ_y(i,j) = y(j) * ( Zi(i,j+1) - Zi(i,j  ) ) * dy(j)
            enddo
!$acc end parallel loop
!$acc parallel loop
        do i=1,imax-1
             j=1
            gradZ_x(i,j) = x(i) * ( Zi(i+1,j) - Zi(i,j  ) ) * dx(i+1)
            gradZ_y(i,j) = y(i) * ( Zi(i,j+1) - Zi(i,j  ) ) * dy(j+1)
        enddo
!$acc end parallel loop
!$acc parallel loop collapse(2) private(i, j)
        do i=2,imax
            do j=2,jmax
            gradZ_x(i,j) = x(i) * ( Zi(i,j) - Zi(i-1,j  ) ) * dx(i)
            gradZ_y(i,j) = y(j) * ( Zi(i,j) - Zi(i  ,j-1) ) * dy(j)
            enddo
        enddo
!$acc END PARALLEL loop
!$acc end data
RETURN
END SUBROUTINE bcZ2
!----------------------------------------------------------------
SUBROUTINE bcZ1(Z_n_tau)
    USE comum
    use openacc
    ! $acc routine
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: Z_n_tau,gradZ_x,gradZ_y
!$acc data present (Z_n_tau,x, y,dx,dy)

    !!! boundary condition

!--------contorno inferior e superior --------
!$acc parallel loop
    do i=1,imax
        Z_n_tau(i,1  ) = Z_n_tau(i,2)
        Z_n_tau(i,jmax) = Z_n_tau(i,jmax-1) + 1.d0*(Z_n_tau(i,jmax-1)-Z_n_tau(i,jmax-2))
    enddo
!$acc end parallel loop
!-------contorno esquerdo e direito -------
!$acc parallel loop
    do j=1,jmax
        Z_n_tau(1  ,j) = Z_n_tau(2,j)
        Z_n_tau(imax,j) = Z_n_tau(imax-1,j)
    enddo
!$acc end parallel loop
       i=1
!$acc parallel loop
            do j=2,jmax-1
            gradZ_x(i,j) = x(i) * ( Z_n_tau(i+1,j) - Z_n_tau(i,j  ) ) * dx(i)
            gradZ_y(i,j) = y(j) * ( Z_n_tau(i,j+1) - Z_n_tau(i,j  ) ) * dy(j)
            enddo
!$acc end parallel loop
!$acc parallel loop
        do i=1,imax-1
             j=1
            gradZ_x(i,j) = x(i) * ( Z_n_tau(i+1,j) - Z_n_tau(i,j  ) ) * dx(i+1)
            gradZ_y(i,j) = y(i) * ( Z_n_tau(i,j+1) - Z_n_tau(i,j  ) ) * dy(j+1)
        enddo
!$acc end parallel loop

!$acc parallel loop collapse(2) private(i, j)
        do i=2,imax
            do j=2,jmax
            gradZ_x(i,j) = x(i) * ( Z_n_tau(i,j) - Z_n_tau(i-1,j  ) ) * dx(i)
            gradZ_y(i,j) = y(j) * ( Z_n_tau(i,j) - Z_n_tau(i  ,j-1) ) * dy(j)
            enddo
        enddo
!$acc END PARALLEL loop
!$acc end data
RETURN
END SUBROUTINE bcZ1
!-------------------------------------------------------------------------

SUBROUTINE bcH(H)
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: H

    !!! boundary condition

!--------contorno inferior e superior --------

    do i=1,imax
        H(i,1  ) = H(i,2)
        H(i,jmax) = H(i,jmax-1)
    enddo

!-------contorno esquerdo e direito -------
    
    do j=1,jmax
        H(1  ,j) = H(2,j)
        H(imax,j) = H(imax-1,j)
    enddo
!!$omp parallel do private(i, j) 
    do i=1,imax
        do j=1,jmax
            if (flag(i,j) .eq. C_B.or.flag(i,j) .eq. C_I) then
                H(i,j) = Lf*Tsup/q + 1.d0/ (S + 1.d0)
            endif
        enddo
    enddo
!!$OMP END PARALLEL DO

RETURN
END SUBROUTINE bcH
