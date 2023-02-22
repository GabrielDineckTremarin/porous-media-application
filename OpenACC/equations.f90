subroutine RESU(um,vm,p,RU)
    use comum
    !$acc routine(alpha)   
     use openacc
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: P
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: q_art, artMAX
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
  !$omp parallel do private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dudxdx,dxdvdy,q_art,u_W,u_WW,u_E,u_EE,u_N,u_NN,u_S,u_SS,v_P,u_P&
  !$OMP&,ap,aw,aww,ae,aee,an,ann,as,ass)
DO i=3,imax-1
    DO j=3,jmax-2

        fn = 0.5d0 * ( vm(i  ,j+1) + vm(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        fs = 0.5d0 * ( vm(i  ,j  ) + vm(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        fe = 0.5d0 * ( um(i+1,j  ) + um(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        fw = 0.5d0 * ( um(i  ,j  ) + um(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )
        df = fe - fw + fn - fs
        Dn = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        Ds = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        De = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        Dw = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))
        aw = Dw + 0.75d0  * alpha(fw) * fw &
                + 0.125d0 * alpha(fe) * fe &
                + 0.375d0 * ( 1.0d0 - alpha(fw) ) * fw
        ae = De - 0.375d0* alpha(fe) * fe &
                - 0.75d0  * ( 1.0d0 - alpha(fe) ) * fe &
                - 0.125d0 * ( 1.0d0 - alpha(fw) ) * fw             
        as = Ds + 0.75d0  * alpha(fs) * fs &
                + 0.125d0 * alpha(fn) * fn &
                + 0.375d0 * ( 1.0d0 - alpha(fs) ) * fs
        an = Dn - 0.375d0* alpha(fn) * fn &
                - 0.75d0  * ( 1.0d0 - alpha(fn) ) * fn & 
                - 0.125d0 * ( 1.0d0 - alpha(fs) ) * fs
        aww = -0.125d0 *           alpha(fw)   * fw
        aee =  0.125d0 * ( 1.0d0 - alpha(fe) ) * fe
        ass = -0.125d0 *           alpha(fs)   * fs
        ann =  0.125d0 * ( 1.0d0 - alpha(fn) ) * fn
        ap = aw + ae + as + an + aww + aee + ass + ann + df
        u_W  = um(i-1,j  )
        u_WW = um(i-2,j  )
        u_E  = um(i+1,j  )
        u_EE = um(i+2,j  )
        u_S  = um(i  ,j-1)
        u_SS = um(i  ,j-2)
        u_N  = um(i  ,j+1)
        u_NN = um(i  ,j+2)        
        u_P = um(i  ,j  )
        v_P = vm(i  ,j  )
            dudxdx =  areau_e(j) * ( u_E - u_P ) / (xm(i+1)-xm(i  )) &
                     -areau_w(j) * ( u_P - u_W ) / (xm(i  )-xm(i-1)) 
            
            dxdvdy =  areau_e(j) * (vm(i  ,j+1) - vm(i  ,j  )) / (ym(j+1)-ym(j)) &
                     -areau_w(j) * (vm(i-1,j+1) - vm(i-1,j  )) / (ym(j+1)-ym(j))
            
            artDivU(i,j) = - b_art *(dudxdx + dxdvdy)   
            
            !bulk artificial viscosity term from Ramshaw(1990)
            q_art =epsilon1(i,j)*  ( p(i,j)-p(i-1,j) )/ (x(i)-x(i-1)) + artDivU(i,j)

        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * u_P &
                                  + aww * u_WW + aw * u_W &
                                  + aee * u_EE + ae * u_E &
                                  + ass * u_SS + as * u_S &
                                  + ann * u_NN + an * u_N) &
                                  - q_art&
               -epsilon1(i,j)*(  u_P/(Re*Darcy_number) + &
    CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* u_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)

    ENDDO
ENDDO    
  !$omp end parallel do

!  $omp parallel do private(i) 
DO i=2,imax
    CALL upwind_U(um,vm,p,RU,i,2)
    CALL upwind_U(um,vm,p,RU,i,jmax-1)
ENDDO
!  $omp end parallel do

!$omp parallel do private(j) 
DO j=2,jmax-1
    CALL upwind_U(um,vm,p,RU,2,j)
    CALL upwind_U(um,vm,p,RU,imax,j)
ENDDO
!$omp end parallel do

RETURN
END SUBROUTINE RESU

!#######################################################################################
subroutine RESU1(um_tau,vm_tau,p,RU, i,j,areau_n, areau_s, areau_e, areau_w, epsilon1, Re, x, y, xm,ym, artDivU, b_art, liga_poros, var) !um,vm,p,RU)
    use comum   
    use openacc
    !$acc routine 
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: P
    REAL(8), DIMENSION(45) :: var
    !real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    !REAL(8) :: Dn, Ds, De, Dw, q_art
    !real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    !real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: artDivU(imax,jmax)
    !real(8) :: afw, afe, afn, afs 
    !real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
    REAL(8)            :: areau_n(2:imax)            !area n de u
    REAL(8)            :: areau_s(2:imax)            !area s de u
    REAL(8)            :: areau_e(1:jmax)            !area e de u
    REAL(8)            :: areau_w(1:jmax)
    REAL(8)            :: epsilon1(imax,jmax)
    REAL(8)            :: Re
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: b_art         !coeficiente da dissipação artificial
    REAL(8)            :: liga_poros(imax,jmax) 

!$acc data present(um_tau,vm_tau,P,RU,areau_n,areau_s,areau_e,areau_w,epsilon1,x, y, xm,ym,artDivU,liga_poros, var)
var = 0.d0
        var(3) = 0.5d0 * ( vm_tau(i  ,j+1) + vm_tau(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        var(2) = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        var(1) = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        var(40) = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )
        var(41) = var(1) - var(40) + var(3) - var(2)
        var(13) = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        var(14) = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        var(15) = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        var(16) = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))
        
    if(var(40).GT.0.0d0) then
       var(42) = 1.0d0
    elseif(var(40).LT.0.0d0) then
       var(42) = 0.0d0
    endif

    if(var(1).GT.0.0d0) then
       var(43) = 1.0d0
    elseif(var(1).LT.0.0d0) then
       var(43) = 0.0d0
    endif

    if(var(3).GT.0.0d0) then
       var(44) = 1.0d0
    elseif(var(3).LT.0.0d0) then
       var(44) = 0.0d0
    endif

    if(var(2).GT.0.0d0) then
       var(45) = 1.0d0
    elseif(var(2).LT.0.0d0) then
       var(45) = 0.0d0
    endif

        var(4) = var(16) + 0.75d0  * var(42) * var(40) + 0.125d0 * var(43) *  var(1) &
               + 0.375d0 * ( 1.0d0 - var(42)) * var(40)
        var(6) = var(15) - 0.375d0 * var(43) * var(1) &
               - 0.75d0  * ( 1.0d0 - var(43) ) * var(1) &
              - 0.125d0 * ( 1.0d0 - var(42) ) * var(40)
        var(8) = var(14) + 0.75d0  * var(45) * var(2) &
               + 0.125d0 * var(44) * var(3) &
               + 0.375d0 * ( 1.0d0 - var(45) ) * var(2)
        var(10) = var(13) - 0.375d0* var(44) * var(3) &
               - 0.75d0  * ( 1.0d0 - var(44) ) * var(3)&
               - 0.125d0 * ( 1.0d0 - var(45) ) * var(2)
        var(5) = -0.125d0 * var(42) * var(40)
        var(7) =  0.125d0 * ( 1.0d0 - var(43) ) * var(1)
        var(9) = -0.125d0 * var(45) * var(2)
        var(11) =  0.125d0 * ( 1.0d0 - var(44) ) * var(3)

        var(12) = var(4) + var(6) + var(8) + var(10) + var(5) + var(7) + var(9) + var(11) + var(41)

        var(18) = um_tau(i-1,j  )
        var(19) = um_tau(i-2,j  )
        var(20) = um_tau(i+1,j  )
        var(21) = um_tau(i+2,j  )
        var(22) = um_tau(i  ,j-1)
        var(23) = um_tau(i  ,j-2)
        var(24) = um_tau(i  ,j+1)
        var(25) = um_tau(i  ,j+2)
        var(26) = um_tau(i  ,j  )
        var(35) = vm_tau(i  ,j  )
            var(36) =  areau_e(j) * (var(20) - var(26) ) / (xm(i+1)-xm(i  )) &
                     -areau_w(j) * ( var(26) - var(18)) / (xm(i  )-xm(i-1))
            var(38) =  areau_e(j) * (vm_tau(i  ,j+1) - vm_tau(i  ,j  )) / (ym(j+1)-ym(j)) &
                     -areau_w(j) * (vm_tau(i-1,j+1) - vm_tau(i-1,j  )) / (ym(j+1)-ym(j))

            artDivU(i,j) = - b_art *(var(36) + var(38))

            !bulk artificial viscosity term from Ramshaw(1990)
            var(17) =epsilon1(i,j)*  ( p(i,j)-p(i-1,j) )/ (x(i)-x(i-1)) + artDivU(i,j)

        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - var(12) * var(26) &
                                  + var(5) * var(19) + var(4) * var(18) &
                                  + var(7) * var(21) + var(6) * var(20) &
                                  + var(9) * var(23) + var(8) * var(22) &
                                  + var(11) * var(25) + var(10) * var(24)) &
                                  - var(17)&
               -epsilon1(i,j)*(  var(26)/(Re*Darcy_number) + &
    CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* var(26)*((var(26)**2.d0 +var(35)**2.d0)**0.5d0))*liga_poros(i,j)

!    ENDDO
! ENDDO
! $acc end parallel loop
!$acc end data


!  $omp parallel do private(i) 
!DO i=2,imax
!    CALL upwind_U1(um_tau,vm_tau,p,RU,i,2)
!    CALL upwind_U1(um_tau,vm_tau,p,RU,i,jmax-1)
!ENDDO
!  $omp end parallel do

! $omp parallel do private(j) 
!DO j=2,jmax-1
!    CALL upwind_U1(um_tau,vm_tau,p,RU,2,j)
!    CALL upwind_U1(um_tau,vm_tau,p,RU,imax,j)
!ENDDO
! $omp end parallel do

RETURN
END SUBROUTINE RESU1

!##########################################################################################################
!##########################################################################################################
subroutine RESU2(ui,vm_tau,p,RU, i, j, areau_n, areau_s, areau_e, areau_w, epsilon1, Re, x, y, xm,ym, artDivU, b_art, liga_poros, var)    
    use comum  
    use openacc
    !$acc routine 
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: ui, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: P
    REAL(8), DIMENSION(45) :: var
    !real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    !REAL(8) :: Dn, Ds, De, Dw
    !real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    !real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: artDivU(imax,jmax)
    !real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
    REAL(8)            :: areau_n(2:imax)            !area n de u
    REAL(8)            :: areau_s(2:imax)            !area s de u
    REAL(8)            :: areau_e(1:jmax)            !area e de u
    REAL(8)            :: areau_w(1:jmax)
    REAL(8)            :: epsilon1(imax,jmax)
    REAL(8)            :: Re
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: b_art         !coeficiente da dissipação artificial
    REAL(8)            :: liga_poros(imax,jmax)

!$acc data present(ui,vm_tau,P,RU,areau_n,areau_s,areau_e,areau_w,epsilon1,x, y, xm,ym,artDivU,liga_poros, var)
var = 0.d0

! $acc parallel loop
!!!!!!!!!!!!1 private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dudxdx,dxdvdy,q_art,u_W,u_WW,u_E,u_EE,u_N,u_NN,u_S,u_SS,v_P,u_P,acc&,ap,aw,aww,ae,aee,an,ann,as,ass)
!DO i=3,imax-1
!    DO j=3,jmax-2
        var(3) = 0.5d0 * ( vm_tau(i  ,j+1) + vm_tau(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        var(2) = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        var(1) = 0.5d0 * ( ui(i+1,j  ) + ui(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        var(40) = 0.5d0 * ( ui(i  ,j  ) + ui(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )
        var(41) = var(1) - var(40) + var(3) - var(2)
        var(13) = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        var(14) = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        var(15) = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        var(16) = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))

    if(var(40).GT.0.0d0) then
       var(42) = 1.0d0
    elseif(var(40).LT.0.0d0) then
       var(42) = 0.0d0
    endif

    if(var(1).GT.0.0d0) then
       var(43) = 1.0d0
    elseif(var(1).LT.0.0d0) then
       var(43) = 0.0d0
    endif

    if(var(3).GT.0.0d0) then
       var(44) = 1.0d0
    elseif(var(3).LT.0.0d0) then
       var(44) = 0.0d0
    endif

    if(var(2).GT.0.0d0) then
       var(45) = 1.0d0
    elseif(var(2).LT.0.0d0) then
       var(45) = 0.0d0
    endif

        var(4) = var(16) + 0.75d0  * var(42) * var(40) + 0.125d0 * var(43) *  var(1) &
               + 0.375d0 * ( 1.0d0 - var(42)) * var(40)
        var(6) = var(15) - 0.375d0 * var(43) * var(1) &
               - 0.75d0  * ( 1.0d0 - var(43) ) * var(1) &
              - 0.125d0 * ( 1.0d0 - var(42) ) * var(40)
        var(8) = var(14) + 0.75d0  * var(45) * var(2) &
               + 0.125d0 * var(44) * var(3) &
               + 0.375d0 * ( 1.0d0 - var(45) ) * var(2)
        var(10) = var(13) - 0.375d0* var(44) * var(3) &
               - 0.75d0  * ( 1.0d0 - var(44) ) * var(3)&
               - 0.125d0 * ( 1.0d0 - var(45) ) * var(2)
        var(5) = -0.125d0 * var(42) * var(40)
        var(7) =  0.125d0 * ( 1.0d0 - var(43) ) * var(1)
        var(9) = -0.125d0 * var(45) * var(2)
        var(11) =  0.125d0 * ( 1.0d0 - var(44) ) * var(3)

        var(12) = var(4) + var(6) + var(8) + var(10) + var(5) + var(7) + var(9) + var(11) + var(41)

        var(18) = ui(i-1,j  )
        var(19) = ui(i-2,j  )
        var(20) = ui(i+1,j  )
        var(21) = ui(i+2,j  )
        var(22) = ui(i  ,j-1)
        var(23) = ui(i  ,j-2)
        var(24) = ui(i  ,j+1)
        var(25) = ui(i  ,j+2)
        var(26) = ui(i  ,j  )
        var(35) = vm_tau(i  ,j  )
            var(36) =  areau_e(j) * (var(20) - var(26) ) / (xm(i+1)-xm(i  )) &
                     -areau_w(j) * ( var(26) - var(18)) / (xm(i  )-xm(i-1))
            var(38) =  areau_e(j) * (vm_tau(i  ,j+1) - vm_tau(i  ,j  )) / (ym(j+1)-ym(j)) &
                     -areau_w(j) * (vm_tau(i-1,j+1) - vm_tau(i-1,j  )) / (ym(j+1)-ym(j))

            artDivU(i,j) = - b_art *(var(36) + var(38))

            !bulk artificial viscosity term from Ramshaw(1990)
            var(17) =epsilon1(i,j)*  ( p(i,j)-p(i-1,j) )/ (x(i)-x(i-1)) + artDivU(i,j)

        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - var(12) * var(26) &
                                  + var(5) * var(19) + var(4) * var(18) &
                                  + var(7) * var(21) + var(6) * var(20) &
                                  + var(9) * var(23) + var(8) * var(22) &
                                  + var(11) * var(25) + var(10) * var(24)) &
                                  - var(17)&
               -epsilon1(i,j)*(  var(26)/(Re*Darcy_number) + &
    CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* var(26)*((var(26)**2.d0 +var(35)**2.d0)**0.5d0))*liga_poros(i,j)
     
!    ENDDO
!ENDDO

!$acc end data

!DO i=2,imax
!    CALL upwind_U2(ui,vm_tau,p,RU,i,2)
!    CALL upwind_U2(ui,vm_tau,p,RU,i,jmax-1)
!ENDDO
!DO j=2,jmax-1
!    CALL upwind_U2(ui,vm_tau,p,RU,2,j)
!    CALL upwind_U2(ui,vm_tau,p,RU,imax,j)
!ENDDO
RETURN
END SUBROUTINE RESU2

!#######################################################################################
!#######################################################################################

subroutine upwind_U(um,vm,p,RU,i,j)
    use comum
    use openacc
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: P
 !   REAL(8), DIMENSION(1:imax+1,1:jmax+1) ::epsilon1,Darcy
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: q_art, artMAX
    !real(8) :: q_art, artDivU(imax,jmax), artDivV(imax,jmax), artMAX
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
        !!!!!!!!!!!!!!!!!!!!!compute x-direction velocity component un!!!!!!!!!!!!!!!!!!!!!
        fn =0.5d0 * ( vm(i  ,j+1) + vm(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        fs = 0.5d0 * ( vm(i  ,j  ) + vm(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        fe = 0.5d0 * ( um(i+1,j  ) + um(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        fw = 0.5d0 * ( um(i  ,j  ) + um(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )
        df = fe - fw + fn - fs
        Dn = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        Ds = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        De = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        Dw = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))
        aw = Dw + MAX(fw , 0.0d0)
        as = Ds + MAX(fs , 0.0d0)
        ae = De + MAX(0.0d0 , -fe)
        an = Dn + MAX(0.0d0 , -fn)
        ap = aw + ae + as + an + df
        u_W = um(i-1,j  )
        u_E = um(i+1,j  )
        u_S = um(i  ,j-1)
        u_N = um(i  ,j+1)
        u_P = um(i  ,j  )
        v_P = vm(i  ,j  )
            dudxdx =  areau_e(j) * ( u_E - u_P ) / (xm(i+1)-xm(i  )) &
                     -areau_w(j) * ( u_P - u_W ) / (xm(i  )-xm(i-1)) 
            
            dxdvdy =  areau_e(j) * (vm(i  ,j+1) - vm(i  ,j  )) / (ym(j+1)-ym(j)) &
                     -areau_w(j) * (vm(i-1,j+1) - vm(i-1,j  )) / (ym(j+1)-ym(j))
!            
!            !bulk artificial viscosity term from Ramshaw(1990)
            q_art = epsilon1(i,j)* ( p(i,j)-p(i-1,j) )/ (x(i)-x(i-1)) - b_art * (dudxdx + dxdvdy)
        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * u_P &
                                  + aw * u_W + ae * u_E  &
                                  + as * u_S + an * u_N) &
                                  - q_art &
               - epsilon1(i,j)* (  u_P/(Re*Darcy_number) +&
        CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* u_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
RETURN
END SUBROUTINE upwind_U
!#######################################################################################

subroutine upwind_U1(um_tau,vm_tau,p,RU,i,j, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
    use comum
    use openacc
    !$acc routine
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: P
    REAL (8), DIMENSION(45) :: var
    !real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    !REAL(8) :: Dn, Ds, De, Dw
    !real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    !real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    !real(8) :: q_art, artMAX
    !real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
    REAL(8)            :: areau_n(2:imax)            !area n de u
    REAL(8)            :: areau_s(2:imax)            !area s de u
    REAL(8)            :: areau_e(1:jmax)            !area e de u
    REAL(8)            :: areau_w(1:jmax)
    REAL(8)            :: epsilon1(imax,jmax)
    REAL(8)            :: Re
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: b_art         !coeficiente da dissipação artificial
    REAL(8)            :: liga_poros(imax,jmax)

!$acc data present (um_tau, vm_tau, P, RU, areau_n,areau_s,areau_e,areau_w,epsilon1,x, y,xm,ym, liga_poros, var)
        var(3) = 0.5d0 * ( vm_tau(i  ,j+1) + vm_tau(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        var(2) = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        var(1) = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        var(40) = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )
        var(41) = var(1) - var(40) + var(3) - var(2)
        var(13) = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        var(14) = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        var(15) = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        var(16) = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))
        var(4) = var(16) + MAX(var(40), 0.0d0)
        var(8) = var(14) + MAX(var(2), 0.0d0)
        var(6) = var(15) + MAX(0.0d0 , - var(1))
        var(10) = var(13) + MAX(0.0d0 , - var(3))
        var(12) = var(4) + var(6) + var(8) + var(10) + var(41)
        var(18) = um_tau(i-1,j  )
        var(20) = um_tau(i+1,j  )
        var(22) = um_tau(i  ,j-1)
        var(24) = um_tau(i  ,j+1)
        var(26) = um_tau(i  ,j  )
        var(35) = vm_tau(i  ,j  )
            var(36) =  areau_e(j) * ( var(20) - var(26) ) / (xm(i+1)-xm(i  )) &
                     -areau_w(j) * ( var(26) - var(18)  ) / (xm(i  )-xm(i-1))

            var(38) =  areau_e(j) * (vm_tau(i  ,j+1) - vm_tau(i  ,j  )) / (ym(j+1)-ym(j)) &
                     -areau_w(j) * (vm_tau(i-1,j+1) - vm_tau(i-1,j  )) / (ym(j+1)-ym(j))
!            
!            !bulk artificial viscosity term from Ramshaw(1990)
            var(17) = epsilon1(i,j)* ( p(i,j)-p(i-1,j) )/ (x(i)-x(i-1)) - b_art * (var(36) + var(38))
	    RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - var(12) * var(26) &
                                  + var(4) * var(18)  + var(6) * var(20)  &
                                  + var(8) * var(22) + var(10) * var(24)) &
                                  - var(17) &
               - epsilon1(i,j)* (  var(26)/(Re*Darcy_number) +&
        CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* var(26)*((var(26)**2.d0 +var(35)**2.d0)**0.5d0))*liga_poros(i,j)
!$acc end data
RETURN
END SUBROUTINE upwind_U1
!##################################################################################################################
!##################################################################################################################

subroutine upwind_U2(ui,vm_tau,p,RU,i,j, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
    use comum
    use openacc
!$acc routine
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: ui, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: P
    REAL(8), DIMENSION(45) :: var
    !real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    !REAL(8) :: Dn, Ds, De, Dw
    !real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    !real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha !, q_art, artMAX, dudxdx, dvdydy, dxdvdy, dydudx
    REAL(8)            :: areau_n(2:imax)            !area n de u
    REAL(8)            :: areau_s(2:imax)            !area s de u
    REAL(8)            :: areau_e(1:jmax)            !area e de u
    REAL(8)            :: areau_w(1:jmax)
    REAL(8)            :: epsilon1(imax,jmax)
    REAL(8)            :: Re
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: b_art         !coeficiente da dissipação artificial
    REAL(8)            :: liga_poros(imax,jmax)    
    
!$acc data present (ui, vm_tau, P, RU,areau_n,areau_s,areau_e,areau_w,epsilon1(imax,jmax),x, y,xm,ym, liga_poros, var)  
 !!!!!!!!!!!!!!!!!!!!!compute x-direction velocity component un!!!!!!!!!!!!!!!!!!!!!

        var(3) = 0.5d0 * ( vm_tau(i  ,j+1) + vm_tau(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        var(2) = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        var(1) = 0.5d0 * ( ui(i+1,j  ) + ui(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        var(40) = 0.5d0 * ( ui(i  ,j  ) + ui(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )
        var(41) = var(1) - var(40) + var(3) - var(2)
        var(13) = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        var(14) = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        var(15) = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        var(16) = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))
        var(4) = var(16) + MAX(var(40), 0.0d0)
        var(8) = var(14) + MAX(var(2), 0.0d0)
        var(6) = var(15) + MAX(0.0d0 , - var(1))
        var(10) = var(13) + MAX(0.0d0 , - var(3))
        var(12) = var(4) + var(6) + var(8) + var(10) + var(41)
        var(18) = ui(i-1,j  )
        var(20) = ui(i+1,j  )
        var(22) = ui(i  ,j-1)
        var(24) = ui(i  ,j+1)
        var(26) = ui(i  ,j  )
        var(35) = vm_tau(i  ,j  )
            var(36) =  areau_e(j) * ( var(20) - var(26) ) / (xm(i+1)-xm(i  )) &
                     -areau_w(j) * ( var(26) - var(18)  ) / (xm(i  )-xm(i-1))

            var(38) =  areau_e(j) * (vm_tau(i  ,j+1) - vm_tau(i  ,j  )) / (ym(j+1)-ym(j)) &
                     -areau_w(j) * (vm_tau(i-1,j+1) - vm_tau(i-1,j  )) / (ym(j+1)-ym(j))
!            
!            !bulk artificial viscosity term from Ramshaw(1990)
            var(17) = epsilon1(i,j)* ( p(i,j)-p(i-1,j) )/ (x(i)-x(i-1)) - b_art * (var(36) + var(38))
	    RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - var(12) * var(26) &
                                  + var(4) * var(18)  + var(6) * var(20)  &
                                  + var(8) * var(22) + var(10) * var(24)) &
                                  - var(17) &
               - epsilon1(i,j)* (  var(26)/(Re*Darcy_number) +&
        CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* var(26)*((var(26)**2.d0 +var(35)**2.d0)**0.5d0))*liga_poros(i,j)
        
!$acc end data
RETURN
END SUBROUTINE upwind_U2

!#######################################################################################
!#######################################################################################

subroutine solve_U(um,vm,um_n,um_tau,vm_tau,um_n_tau,p,residual_u, var,res_u,RU,UI,areau_e,areau_w,areau_n,areau_s,ym,xm,Re,x,y,liga_poros,dtau,epsilon1,artDivU)
    use comum
    use openacc
    !acc routine
    !$acc routine (RESU1)
    !$acc routine (RESU2)
    !$acc routine (upwind_U1)
    !$acc routine (upwind_U2)
    !$acc routine (bcUV1)
    !$acc routine (bcUV22)
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um,um_n, res_u
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax  ) :: p
    REAL(8), DIMENSION(45) :: var
    REAL(8) :: residual_u
    REAL(8), DIMENSION(1:imax+1,1:jmax) :: RU, ui
    REAL(8)            :: artDivU(imax,jmax)
    REAL(8)            :: areau_n(2:imax)            !area n de u
    REAL(8)            :: areau_s(2:imax)            !area s de u
    REAL(8)            :: areau_e(1:jmax)            !area e de u
    REAL(8)            :: areau_w(1:jmax)            !area w de u
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: Re,dtau
    REAL(8)            :: liga_poros(imax,jmax)
    REAL(8)            :: epsilon1(imax,jmax)
  
    RU = 0.d0
    res_u = 0.d0
    ui = 0.d0
    var = 0.d0
! $acc enter data copyin(um, um_n, vm, um_tau, um_n_tau, vm_tau, p, residual_u, RU, res_u, ui,&
! $acc&areau_n, areau_s, areau_e, areau_w,epsilon1, x, y, xm,ym, b_art, liga_poros, Re, artDivU, v_i, var)

!$acc parallel loop collapse(2) private(var)
DO i=3,imax-1
    DO j=3,jmax-2
       CALL RESU1(um_tau,vm_tau,p,RU,i,j,areau_n, areau_s, areau_e, areau_w, epsilon1, Re, x, y, xm,ym, artDivU, b_art, liga_poros, var)
    ENDDO
ENDDO 
!$acc end parallel loop    

!$acc parallel loop private(var, i)
    DO i=2,imax
       CALL upwind_U1(um_tau,vm_tau,p,RU,i,2, areau_n, areau_s, areau_e, areau_w,epsilon1,Re,x,y,xm,ym,b_art,liga_poros, var)
       CALL upwind_U1(um_tau,vm_tau,p,RU,i,jmax-1, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
    ENDDO
!$acc end parallel loop
!$acc parallel loop private(j, var)
    DO j=2,jmax-1
      CALL upwind_U1(um_tau,vm_tau,p,RU,2,j, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
      CALL upwind_U1(um_tau,vm_tau,p,RU,imax,j, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
    ENDDO
!$acc end parallel loop

!$acc parallel loop collapse(2) private(i,j) 
    DO i=3,imax-1
        DO j=2,jmax-1
            res_u(i,j) =( (um(i,j)-um_tau(i,j)) +  RU(i,j)*dt) * dtau
            ui(i,j) =  ( um_tau(i,j) + res_u(i,j))
        ENDDO
    ENDDO
!$acc end parallel loop

   call bcUV1(ui,vm_tau)

!$acc parallel loop collapse(2) private(var)
DO i=3,imax-1
    DO j=3,jmax-2
          CALL RESU2(ui,vm_tau,p,RU, i, j, areau_n, areau_s, areau_e, areau_w, epsilon1, Re, x, y, xm,ym, artDivU, b_art, liga_poros, var)       
    ENDDO
ENDDO   
!$acc end parallel loop       

!$acc parallel loop private(var)
    DO i=2,imax
       CALL upwind_U2(ui,vm_tau,p,RU,i,2, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
       CALL upwind_U2(ui,vm_tau,p,RU,i,jmax-1, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
    ENDDO
!$acc end parallel loop
!$acc parallel loop private(var)
    DO j=2,jmax-1
       CALL upwind_U2(ui,vm_tau,p,RU,2,j, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
       CALL upwind_U2(ui,vm_tau,p,RU,imax,j, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
    ENDDO
!$acc end parallel loop

    !$acc parallel loop collapse(2) private(i,j) 
    DO i=3,imax-1
        DO j=2,jmax-1
            res_u(i,j) =( (um(i,j)-um_tau(i,j)) +  RU(i,j)*dt) * dtau   
            ui(i,j) = (0.75d0 * um_tau(i,j) + 0.25d0 * &
                      ( ui(i,j) + res_u(i,j) ) )
        ENDDO
    ENDDO
    !$acc end parallel loop

    call bcUV1(ui,vm_tau)

!$acc parallel loop collapse(2) private(var)
DO i=3,imax-1
    DO j=3,jmax-2
	CALL RESU2(ui,vm_tau,p,RU, i, j, areau_n, areau_s, areau_e, areau_w, epsilon1, Re, x, y, xm,ym, artDivU, b_art, liga_poros, var)
    ENDDO
ENDDO    
!$acc end parallel loop
          
!$acc parallel loop private(var)    
    DO i=2,imax
       CALL upwind_U2(ui,vm_tau,p,RU,i,2, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
       CALL upwind_U2(ui,vm_tau,p,RU,i,jmax-1,areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,b_art,liga_poros, var)
    ENDDO
!$acc end parallel loop 
!$acc parallel loop private(var)    
   DO j=2,jmax-1
       CALL upwind_U2(ui,vm_tau,p,RU,2,j, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,artDivU,b_art,liga_poros, var)
       CALL upwind_U2(ui,vm_tau,p,RU,imax,j, areau_n, areau_s, areau_e,areau_w,epsilon1, Re, x, y,xm,ym,artDivU,b_art,liga_poros, var)
    ENDDO
!$acc end parallel loop
    
!$acc parallel loop collapse(2) private(i,j) 
    DO i=3,imax-1
        DO j=2,jmax-1
            res_u(i,j) =( (um(i,j)-um_tau(i,j)) +  RU(i,j)*dt) * dtau
            um_n_tau(i,j) = (1.0d0 / 3.0d0 * um_tau(i,j) + 2.0d0 / 3.0d0 * &
                      ( ui(i,j) + res_u(i,j) ) )  
  	ENDDO
    ENDDO
!$acc end parallel loop

    call bcUV2(um_n_tau, vm_tau)

!$acc parallel
    residual_u =  MAXVAL(ABS(res_u)) 
!$acc end parallel

! $acc exit data copyout(um, um_n, vm, um_tau, um_n_tau, vm_tau, p, residual_u, RU, res_u, ui,&
! $acc&areau_n, areau_s, areau_e, areau_w,epsilon1, x, y, xm,ym, b_art, liga_poros, Re, artDivU, v_i, var)

! $delete(um, um_n, vm, um_tau, um_n_tau, vm_tau, p, residual_u, RU, res_u, ui,&
! $delete&areau_n, areau_s, areau_e, areau_w,epsilon1, x, y, xm,ym, b_art, liga_poros, Re, artDivU, v_i, var)


!    open (550,file='data/resu.dat')!        
!         do i = 1,imax
!             do j = 1,jmax!            
!             write (550,*) X(i) , Y(j) , res_u(i,j)!
!             enddo
!         enddo    
!    close(550)
RETURN
END SUBROUTINE solve_U

!#######################################################################################
!#######################################################################################
!#######################################################################################

subroutine RESV(um,vm,p,T,InvFr2,RV)
    !$acc routine(alpha)
    !$acc routine(upwind_V)
    use comum
    use openacc
    !use omp_lib
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
 !   REAL(8), DIMENSION(1:imax+1,1:jmax+1) :: epsilon1, Darcy
    real(8) :: InvFr2
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: q_art, artMAX
    !real(8) :: q_art, artDivU(imax,jmax), artDivV(imax,jmax), artMAX
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx

DO i=3,imax-2
    DO j=3,jmax-1

        fn =  0.5d0 * ( vm(i  ,j  ) + vm(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        fs = 0.5d0 * ( vm(i  ,j  ) + vm(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        fe = 0.5d0 * ( um(i+1,j  ) + um(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        fw = 0.5d0 * ( um(i  ,j  ) + um(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)

        df = fe - fw + fn - fs

        Dn = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))

        !CALL quick(ap,aw,aww,ae,aee,an,ann,as,ass,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
        !!! init quick
        aw = Dw + 0.75d0  * alpha(fw) * fw &
                + 0.125d0 * alpha(fe) * fe &
                + 0.375d0 * ( 1.0d0 - alpha(fw) ) * fw
        ae = De - 0.375d0* alpha(fe) * fe &
                - 0.75d0  * ( 1.0d0 - alpha(fe) ) * fe &
                - 0.125d0 * ( 1.0d0 - alpha(fw) ) * fw             
        as = Ds + 0.75d0  * alpha(fs) * fs &
                + 0.125d0 * alpha(fn) * fn &
                + 0.375d0 * ( 1.0d0 - alpha(fs) ) * fs
        an = Dn - 0.375d0* alpha(fn) * fn &
                - 0.75d0  * ( 1.0d0 - alpha(fn) ) * fn & 
                - 0.125d0 * ( 1.0d0 - alpha(fs) ) * fs
        aww = -0.125d0 *           alpha(fw)   * fw
        aee =  0.125d0 * ( 1.0d0 - alpha(fe) ) * fe
        ass = -0.125d0 *           alpha(fs)   * fs
        ann =  0.125d0 * ( 1.0d0 - alpha(fn) ) * fn

        ap = aw + ae + as + an + aww + aee + ass + ann + df
        !!! end quick 

        v_W  = vm(i-1,j  )
        v_WW = vm(i-2,j  )
        v_E  = vm(i+1,j  )
        v_EE = vm(i+2,j  )
        v_S  = vm(i  ,j-1)
        v_SS = vm(i  ,j-2)
        v_N  = vm(i  ,j+1)
        v_NN = vm(i  ,j+2)         
        v_P  = vm(i  ,j  )
	u_P  = um(i   ,j  )

            dvdydy =  areav_n(i) * ( v_N - v_P ) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( v_P - v_S ) / (ym(j  )-ym(j-1)) 
            
            dydudx =  areav_n(i) * (um(i+1,j  )-um(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um(i+1,j-1)-um(i,j-1)) / (xm(i+1)-xm(i))

            artDivV(i,j) = - b_art *(dydudx + dvdydy)            

            !bulk artificial viscosity term from Ramshaw(1990)
            q_art = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) + artDivV(i,j)

 
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * (- ap * v_P &
                                  + aww * v_WW + aw * v_W &
                                  + aee * v_EE + ae * v_E &
                                  + ass * v_SS + as * v_S &
                                  + ann * v_NN + an * v_N) &
                                  - q_art &
              + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
               -epsilon1(i,j)*(   v_P/(Re*Darcy_number) + &
     CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
ENDDO    
ENDDO
! $omp end parallel do

!$omp parallel do private (i)
DO i=2,imax-1
    CALL upwind_V(um,vm,p,RV,InvFr2,T,i,2)
    CALL upwind_V(um,vm,p,RV,InvFr2,T,i,jmax)
ENDDO
!$omp end parallel do
!$omp parallel do private (j)
DO j=2,jmax
    CALL upwind_V(um,vm,p,RV,InvFr2,T,2,j)
    CALL upwind_V(um,vm,p,RV,InvFr2,T,imax-1,j)
ENDDO
!$omp end parallel do
RETURN
END SUBROUTINE RESV


!#######################################################################################

subroutine RESV1(um_tau,vm_tau,p,T,InvFr2,RV,i,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,artDivV,b_art,liga_poros,var)
    use comum
    use openacc
    !$acc routine
    !$acc routine(alpha)
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
    REAL(8), DIMENSION(45) :: var
    real(8) :: InvFr2
    real(8) :: artDivV(imax,jmax)
    !REAL(8)            :: afw, afe, afn, afs 
    REAL(8)            :: areav_n(1:imax)            !area n de v
    REAL(8)            :: areav_s(1:imax)            !area s de v
    REAL(8)            :: areav_e(2:jmax)            !area e de v
    REAL(8)            :: areav_w(2:jmax)            !area w de v
    REAL(8)            :: epsilon1(imax,jmax)
    REAL(8)            :: Re
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: b_art         !coeficiente da dissipação artificial
    REAL(8)            :: liga_poros(imax,jmax)

!$acc data present(vi, RV, res_v,um,vm, vm_n, um_tau,vm_tau,vm_n_tau,p,t, v_i, areav_e, areav_w,areav_s,&
!$acc&areav_n, epsilon1, re, ym, x, xm, y, b_art, liga_poros, artDivV, var)
var = 0.d0
!DO i=3,imax-2
    !DO j=3,jmax-1
        var(3)=  0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        var(2) = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        var(1) = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        var(40) = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)
        var(41) = var(1) - var(40) + var(3) - var(2)
        var(13) = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        var(14) = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        var(15) = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        var(16) = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))

            if(var(40).GT.0.0d0) then
               var(42) = 1.0d0
            elseif(var(40).LT.0.0d0) then
               var(42) = 0.0d0
            endif

            if(var(1).GT.0.0d0) then
               var(43) = 1.0d0
            elseif(var(1).LT.0.0d0) then
               var(43) = 0.0d0
            endif

            if(var(3).GT.0.0d0) then
               var(44) = 1.0d0
	    elseif(var(3).LT.0.0d0) then
               var(44) = 0.0d0
            endif

            if(var(2).GT.0.0d0) then
               var(45) = 1.0d0
            elseif(var(2).LT.0.0d0) then
               var(45) = 0.0d0
            endif
        var(4) = var(16) + 0.75d0  * var(42) * var(40) &
                + 0.125d0 * var(43) * var(1) &
                + 0.375d0 * ( 1.0d0 - var(42) ) * var(40)
        var(6) = var(15) - 0.375d0* var(43) * var(1) &
                - 0.75d0  * ( 1.0d0 - var(43) ) * var(1) &
                - 0.125d0 * ( 1.0d0 - var(42) ) * var(40)
        var(8) = var(14) + 0.75d0  * var(45) * var(2) &
                + 0.125d0 * var(44) * var(3) &
                + 0.375d0 * ( 1.0d0 - var(45) ) * var(2)
        var(10) = var(13) - 0.375d0* var(44) * var(3)  &
                - 0.75d0  * ( 1.0d0 - var(44) ) * var(3)  &
                - 0.125d0 * ( 1.0d0 - var(45) ) * var(2)
        var(5) = -0.125d0 *           var(42)   * var(40)
        var(7) =  0.125d0 * ( 1.0d0 - var(43) ) * var(1)
        var(9)  = -0.125d0 *           var(45)   * var(2)
        var(11) =  0.125d0 * ( 1.0d0 - var(44) ) * var(3)
        var(12) = var(4) + var(6) + var(8) + var(10) + var(5) + var(7) + var(9) + var(11) + var(41)
        var(27) = vm_tau(i-1,j  )
        var(28) = vm_tau(i-2,j  )
        var(29) = vm_tau(i+1,j  )
        var(30) = vm_tau(i+2,j  )
        var(31) = vm_tau(i  ,j-1)
        var(32) = vm_tau(i  ,j-2)
        var(33) = vm_tau(i  ,j+1)
        var(34) = vm_tau(i  ,j+2)
        var(35) = vm_tau(i  ,j  )
        var(26) = um_tau(i   ,j  )

            var(37) =  areav_n(i) * ( var(33) - var(35) ) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( var(35) - var(31) ) / (ym(j  )-ym(j-1))

            var(39) =  areav_n(i) * (um_tau(i+1,j  )-um_tau(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um_tau(i+1,j-1)-um_tau(i,j-1)) / (xm(i+1)-xm(i))

            artDivV(i,j) = - b_art *(var(39) + var(37))

            !bulk artificial viscosity term from Ramshaw(1990)
            var(17) = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) + artDivV(i,j)
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * (- var(12) * var(35) &
                                  + var(5) * var(28) + var(4) * var(27) &
                                  + var(7) * var(30) + var(6) * var(29) &
                                  + var(9) * var(32) + var(8) * var(31) &
                                  + var(11) * var(34) + var(10) * var(33)) &
                                  - var(17) &
              + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
               -epsilon1(i,j)*(   var(35)/(Re*Darcy_number) + &
     CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* var(35)*((var(26)**2.d0 +var(35)**2.d0)**0.5d0))*liga_poros(i,j)
!ENDDO
!ENDDO

!$acc end data

RETURN
END SUBROUTINE RESV1

!--------------------------------------------------------------------------------------------------
subroutine RESV2(um_tau,vi,p,T,InvFr2,RV,i,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,artDivV,b_art,liga_poros,var)
    use comum
    use openacc
    !$acc routine
    !$acc routine(alpha)
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    !REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, RV
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vi, RV

    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
    REAL(8), DIMENSION(45) :: var
    real(8) :: InvFr2
    real(8) :: artDivV(imax,jmax)
    REAL(8)            :: areav_n(1:imax)            !area n de v
    REAL(8)            :: areav_s(1:imax)            !area s de v
    REAL(8)            :: areav_e(2:jmax)            !area e de v
    REAL(8)            :: areav_w(2:jmax)            !area w de v
    REAL(8)            :: epsilon1(imax,jmax)
    REAL(8)            :: Re
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: b_art         !coeficiente da dissipação artificial
    REAL(8)            :: liga_poros(imax,jmax)

! $acc data present(vi, RV, res_v,um,vm, vm_n, um_tau,vm_tau,vm_n_tau,p,t, v_i, areav_e, areav_w,areav_s,&
! $acc&areav_n, epsilon1, re, ym, x, xm, y, b_art, liga_poros, artDivV, var)

! $acc parallel private(var, i, j)
!DO i=3,imax-2
!    DO j=3,jmax-1
        var(3)=  0.5d0 * ( vi(i  ,j  ) + vi(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        var(2) = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        var(1) = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        var(40) = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)
        var(41) = var(1) - var(40) + var(3) - var(2)
        var(13) = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        var(14) = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        var(15) = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        var(16) = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))

            if(var(40).GT.0.0d0) then
               var(42) = 1.0d0
            elseif(var(40).LT.0.0d0) then
               var(42) = 0.0d0
            endif

            if(var(1).GT.0.0d0) then
               var(43) = 1.0d0
            elseif(var(1).LT.0.0d0) then
               var(43) = 0.0d0
            endif
            if(var(3).GT.0.0d0) then
               var(44) = 1.0d0
            elseif(var(3).LT.0.0d0) then
               var(44) = 0.0d0
            endif
            if(var(2).GT.0.0d0) then
               var(45) = 1.0d0
            elseif(var(2).LT.0.0d0) then
               var(45) = 0.0d0
            endif
	var(4) = var(16) + 0.75d0  * var(42) * var(40) &
                + 0.125d0 * var(43) * var(1) &
                + 0.375d0 * ( 1.0d0 - var(42) ) * var(40)
        var(6) = var(15) - 0.375d0* var(43) * var(1) &
                - 0.75d0  * ( 1.0d0 - var(43) ) * var(1) &
                - 0.125d0 * ( 1.0d0 - var(42) ) * var(40)
        var(8) = var(14) + 0.75d0  * var(45) * var(2) &
                + 0.125d0 * var(44) * var(3) &
                + 0.375d0 * ( 1.0d0 - var(45) ) * var(2)
        var(10) = var(13) - 0.375d0* var(44) * var(3)  &
                - 0.75d0  * ( 1.0d0 - var(44) ) * var(3)  &
                - 0.125d0 * ( 1.0d0 - var(45) ) * var(2)
        var(5) = -0.125d0 *           var(42)   * var(40)
        var(7) =  0.125d0 * ( 1.0d0 - var(43) ) * var(1)
        var(9)  = -0.125d0 *           var(45)   * var(2)
        var(11) =  0.125d0 * ( 1.0d0 - var(44) ) * var(3)
        var(12) = var(4) + var(6) + var(8) + var(10) + var(5) + var(7) + var(9) + var(11) + var(41)
        var(27) = vi(i-1,j  )
        var(28) = vi(i-2,j  )
        var(29) = vi(i+1,j  )
        var(30) = vi(i+2,j  )
        var(31) = vi(i  ,j-1)
        var(32) = vi(i  ,j-2)
        var(33) = vi(i  ,j+1)
        var(34) = vi(i  ,j+2)
        var(35) = vi(i  ,j  )
        var(26) = um_tau(i   ,j  )

        var(37) =  areav_n(i) * ( var(33) - var(35) ) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( var(35) - var(31) ) / (ym(j  )-ym(j-1))

            var(39) =  areav_n(i) * (um_tau(i+1,j  )-um_tau(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um_tau(i+1,j-1)-um_tau(i,j-1)) / (xm(i+1)-xm(i))

            artDivV(i,j) = - b_art *(var(39) + var(37))

            !bulk artificial viscosity term from Ramshaw(1990)
            var(17) = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) + artDivV(i,j)
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * (- var(12) * var(35) &
                                  + var(5) * var(28) + var(4) * var(27) &
                                  + var(7) * var(30) + var(6) * var(29) &
                                  + var(9) * var(32) + var(8) * var(31) &
                                  + var(11) * var(34) + var(10) * var(33)) &
                                  - var(17) &
              + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
               -epsilon1(i,j)*(   var(35)/(Re*Darcy_number) + &
     CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* var(35)*((var(26)**2.d0 +var(35)**2.d0)**0.5d0))*liga_poros(i,j)

!ENDDO
!ENDDO
! $acc end parallel

! $acc end data

RETURN
END SUBROUTINE RESV2

!---------------------------------------------------------------------------------------------------
subroutine RESV1x(um_tau,vm_tau,p,T,InvFr2,RV)
    !$acc routine(upwind_V)
    use comum
    use openacc
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
    real(8) :: InvFr2
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: q_art, artMAX
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx

DO i=3,imax-2
    DO j=3,jmax-1
        fn =  0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        fs = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        fe = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        fw = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)
        df = fe - fw + fn - fs
        Dn = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))
        aw = Dw + 0.75d0  * alpha(fw) * fw &
                + 0.125d0 * alpha(fe) * fe &
                + 0.375d0 * ( 1.0d0 - alpha(fw) ) * fw
        ae = De - 0.375d0* alpha(fe) * fe &
                - 0.75d0  * ( 1.0d0 - alpha(fe) ) * fe &
                - 0.125d0 * ( 1.0d0 - alpha(fw) ) * fw
        as = Ds + 0.75d0  * alpha(fs) * fs &
                + 0.125d0 * alpha(fn) * fn &
                + 0.375d0 * ( 1.0d0 - alpha(fs) ) * fs
        an = Dn - 0.375d0* alpha(fn) * fn &
                - 0.75d0  * ( 1.0d0 - alpha(fn) ) * fn &
                - 0.125d0 * ( 1.0d0 - alpha(fs) ) * fs
        aww = -0.125d0 *           alpha(fw)   * fw
        aee =  0.125d0 * ( 1.0d0 - alpha(fe) ) * fe
        ass = -0.125d0 *           alpha(fs)   * fs
        ann =  0.125d0 * ( 1.0d0 - alpha(fn) ) * fn
        ap = aw + ae + as + an + aww + aee + ass + ann + df
        v_W  = vm_tau(i-1,j  )
        v_WW = vm_tau(i-2,j  )
        v_E  = vm_tau(i+1,j  )
        v_EE = vm_tau(i+2,j  )
        v_S  = vm_tau(i  ,j-1)
        v_SS = vm_tau(i  ,j-2)
        v_N  = vm_tau(i  ,j+1)
        v_NN = vm_tau(i  ,j+2)
        v_P  = vm_tau(i  ,j  )
        u_P  = um_tau(i   ,j  )

            dvdydy =  areav_n(i) * ( v_N - v_P ) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( v_P - v_S ) / (ym(j  )-ym(j-1))

            dydudx =  areav_n(i) * (um_tau(i+1,j  )-um_tau(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um_tau(i+1,j-1)-um_tau(i,j-1)) / (xm(i+1)-xm(i))

            artDivV(i,j) = - b_art *(dydudx + dvdydy)

            !bulk artificial viscosity term from Ramshaw(1990)
            q_art = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) + artDivV(i,j)
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * (- ap * v_P &
                                  + aww * v_WW + aw * v_W &
                                  + aee * v_EE + ae * v_E &
                                  + ass * v_SS + as * v_S &
                                  + ann * v_NN + an * v_N) &
                                  - q_art &
              + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
               -epsilon1(i,j)*(   v_P/(Re*Darcy_number) + &
     CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
ENDDO
ENDDO
! $omp end parallel do
! $omp parallel do private (i)
!DO i=2,imax-1
!    CALL upwind_V1(um_tau,vm_tau,p,RV,InvFr2,T,i,2)
!    CALL upwind_V1(um_tau,vm_tau,p,RV,InvFr2,T,i,jmax)
!ENDDO
! $omp end parallel do
! $omp parallel do private (j)
!DO j=2,jmax
!    CALL upwind_V1(um_tau,vm_tau,p,RV,InvFr2,T,2,j)
!    CALL upwind_V1(um_tau,vm_tau,p,RV,InvFr2,T,imax-1,j)
!ENDDO
! $omp end parallel do

RETURN
END SUBROUTINE RESV1x

!#######################################################################################
!#######################################################################################
subroutine RESV2x(um_tau,vi,p,T,InvFr2,RV)
    !$acc routine(upwind_V)
    use comum
    use openacc
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vi, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
    real(8) :: InvFr2
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: q_art, artMAX
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
DO i=3,imax-2
    DO j=3,jmax-1
        fn =  0.5d0 * ( vi(i  ,j  ) + vi(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        fs = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        fe = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        fw = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)
        df = fe - fw + fn - fs
        Dn = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))
        aw = Dw + 0.75d0  * alpha(fw) * fw &
                + 0.125d0 * alpha(fe) * fe &
                + 0.375d0 * ( 1.0d0 - alpha(fw) ) * fw
        ae = De - 0.375d0* alpha(fe) * fe &
                - 0.75d0  * ( 1.0d0 - alpha(fe) ) * fe &
                - 0.125d0 * ( 1.0d0 - alpha(fw) ) * fw
        as = Ds + 0.75d0  * alpha(fs) * fs &
                + 0.125d0 * alpha(fn) * fn &
                + 0.375d0 * ( 1.0d0 - alpha(fs) ) * fs
        an = Dn - 0.375d0* alpha(fn) * fn &
                - 0.75d0  * ( 1.0d0 - alpha(fn) ) * fn &
                - 0.125d0 * ( 1.0d0 - alpha(fs) ) * fs
        aww = -0.125d0 *           alpha(fw)   * fw
        aee =  0.125d0 * ( 1.0d0 - alpha(fe) ) * fe
        ass = -0.125d0 *           alpha(fs)   * fs
        ann =  0.125d0 * ( 1.0d0 - alpha(fn) ) * fn
        ap = aw + ae + as + an + aww + aee + ass + ann + df
        v_W  = vi(i-1,j  )
        v_WW = vi(i-2,j  )
        v_E  = vi(i+1,j  )
        v_EE = vi(i+2,j  )
        v_S  = vi(i  ,j-1)
        v_SS = vi(i  ,j-2)
        v_N  = vi(i  ,j+1)
        v_NN = vi(i  ,j+2)
        v_P  = vi(i  ,j  )
        u_P  = um_tau(i   ,j  )

            dvdydy =  areav_n(i) * ( v_N - v_P ) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( v_P - v_S ) / (ym(j  )-ym(j-1))

            dydudx =  areav_n(i) * (um_tau(i+1,j  )-um_tau(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um_tau(i+1,j-1)-um_tau(i,j-1)) / (xm(i+1)-xm(i))

            artDivV(i,j) = - b_art *(dydudx + dvdydy)

            !bulk artificial viscosity term from Ramshaw(1990)
            q_art = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) + artDivV(i,j)
 RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * (- ap * v_P &
                                  + aww * v_WW + aw * v_W &
                                  + aee * v_EE + ae * v_E &
                                  + ass * v_SS + as * v_S &
                                  + ann * v_NN + an * v_N) &
                                  - q_art &
              + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
               -epsilon1(i,j)*(   v_P/(Re*Darcy_number) + &
     CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
ENDDO
ENDDO
! $omp end parallel do

! $omp parallel do private (i)
!DO i=2,imax-1
!    CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,i,2)
!    CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,i,jmax)
!ENDDO
! $omp end parallel do
! $omp parallel do private (j)
!DO j=2,jmax
!    CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,2,j)
!    CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,imax-1,j)
!ENDDO
! $omp end parallel do

RETURN
END SUBROUTINE RESV2x

!#######################################################################################
!#######################################################################################

subroutine upwind_V2(um_tau,vi,p,RV,InvFr2,T,i,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
    use openacc
    use comum
    !$acc routine
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vi, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
    REAL(8), DIMENSION(45) :: var
    REAL(8) :: InvFr2
    REAL(8)            :: areav_n(1:imax)            !area n de v
    REAL(8)            :: areav_s(1:imax)            !area s de v
    REAL(8)            :: areav_e(2:jmax)            !area e de v
    REAL(8)            :: areav_w(2:jmax)            !area w de v   
    REAL(8)            :: epsilon1(imax,jmax)
    REAL(8)            :: Re
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: b_art         !coeficiente da dissipação artificial
    REAL(8)            :: liga_poros(imax,jmax)
        var(3) = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        var(2) = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        var(1) = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        var(40)= 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)
        var(41) = var(1) - var(40) + var(3) - var(2)
        var(13) = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        var(14) = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        var(15) = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        var(16) = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))
        var(4)  = var(16) + MAX(var(40) , 0.0d0)
        var(8)  = var(14) + MAX(var(2) , 0.0d0)
        var(6)  = var(15) + MAX(0.0d0 , - var(1))
        var(10) = var(13) + MAX(0.0d0 , - var(3))
        var(12) = var(4) + var(6) + var(8) + var(10) + var(41)
        var(27) = vi(i-1,j  )
        var(29) = vi(i+1,j  )
        var(31) = vi(i  ,j-1)
        var(33) = vi(i  ,j+1)
        var(35) = vi(i  ,j  )
        var(26) = um_tau(i  ,j  )     !!!!!POROUS NEW
        var(37) =  areav_n(i) * ( var(33) - var(35)) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( var(35) - var(31)) / (ym(j  )-ym(j-1))

        var(39) =  areav_n(i) * (um_tau(i+1,j  )-um_tau(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um_tau(i+1,j-1)-um_tau(i,j-1)) / (xm(i+1)-xm(i))
            !bulk artificial viscosity term from Ramshaw(1990)
        var(17) = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) - b_art * (var(39) + var(37))
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - var(12) * var(35) &
                                  + var(4) * var(27) + var(6) * var(29)  &
                                  + var(8) * var(31) + var(10) * var(33))  &
                                  - var(17) &
              + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
               -epsilon1(i,j)*  ( var(35)/(Re*Darcy_number) + &
              CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* var(35)*((var(26)**2.d0 +var(35)**2.d0)**0.5d0))*liga_poros(i,j)
              
return
end subroutine upwind_V2
!-------------------------------------------------------------------------------------------

subroutine upwind_V(um,vm,p,RV,InvFr2,T,i,j)
    use openacc
    use comum
    implicit none
    integer :: i, j
    
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
!    REAL(8), DIMENSION(1:imax+1  ,1:jmax+1) :: epsilon1,Darcy
    real(8) :: InvFr2
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: q_art, artMAX
    !real(8) :: q_art, artDivU(imax,jmax), artDivV(imax,jmax), artMAX
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx


        fn = 0.5d0 * ( vm(i  ,j  ) + vm(i  ,j+1) ) * areav_n(i) / epsilon1(i,j) 
        fs = 0.5d0 * ( vm(i  ,j  ) + vm(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)  
        fe = 0.5d0 * ( um(i+1,j  ) + um(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)  
        fw = 0.5d0 * ( um(i  ,j  ) + um(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)  
        df = fe - fw + fn - fs
        Dn = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))
        !WCALL upwind(ap,aw,ae,an,as,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
!upwind 
        aw = Dw + MAX(fw , 0.0d0)
        as = Ds + MAX(fs , 0.0d0)
        ae = De + MAX(0.0d0 , -fe)
        an = Dn + MAX(0.0d0 , -fn)
        ap = aw + ae + as + an + df      
!upwind
	v_W = vm(i-1,j  )
        v_E = vm(i+1,j  )
        v_S = vm(i  ,j-1)
        v_N = vm(i  ,j+1)
        v_P = vm(i  ,j  )
        u_P = um(i  ,j  )     !!!!!POROUS NEW

            dvdydy =  areav_n(i) * ( v_N - v_P ) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( v_P - v_S ) / (ym(j  )-ym(j-1)) 
            
            dydudx =  areav_n(i) * (um(i+1,j  )-um(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um(i+1,j-1)-um(i,j-1)) / (xm(i+1)-xm(i))
            
            !bulk artificial viscosity term from Ramshaw(1990)
            q_art = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) - b_art * (dydudx + dvdydy)
   

 
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * v_P &
                                  + aw * v_W + ae * v_E  &
                                  + as * v_S + an * v_N)  &
                                  - q_art &
              + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
               -epsilon1(i,j)*  ( v_P/(Re*Darcy_number) + &
              CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)

return
end subroutine upwind_V

!#######################################################################################
subroutine upwind_V1(um_tau,vm_tau,p,RV,InvFr2,T,i,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
    use openacc
    use comum
    !$acc routine
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
    REAL(8), DIMENSION(45) :: var
    REAL(8) :: InvFr2
    REAL(8)            :: areav_n(1:imax)            !area n de v
    REAL(8)            :: areav_s(1:imax)            !area s de v
    REAL(8)            :: areav_e(2:jmax)            !area e de v
    REAL(8)            :: areav_w(2:jmax)            !area w de v   
    REAL(8)            :: epsilon1(imax,jmax)
    REAL(8)            :: Re
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: b_art         !coeficiente da dissipação artificial
    REAL(8)            :: liga_poros(imax,jmax)
!$acc data present( RV, um_tau,vm_tau,p,t,areav_e, areav_w,areav_s,areav_n, epsilon1,re,ym,x,xm,y,b_art, liga_poros,var, residual_v)

        var(3) = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        var(2) = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        var(1) = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        var(40)= 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)
        var(41) = var(1) - var(40) + var(3) - var(2)
        var(13) = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        var(14) = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        var(15) = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        var(16) = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))
        var(4)  = var(16) + MAX(var(40) , 0.0d0)
        var(8)  = var(14) + MAX(var(2) , 0.0d0)
        var(6)  = var(15) + MAX(0.0d0 , - var(1))
        var(10) = var(13) + MAX(0.0d0 , - var(3))
        var(12) = var(4) + var(6) + var(8) + var(10) + var(41)
        var(27) = vm_tau(i-1,j  )
        var(29) = vm_tau(i+1,j  )
        var(31) = vm_tau(i  ,j-1)
        var(33) = vm_tau(i  ,j+1)
        var(35) = vm_tau(i  ,j  )
        var(26) = um_tau(i  ,j  )     !!!!!POROUS NEW
        var(37) =  areav_n(i) * ( var(33) - var(35)) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( var(35) - var(31)) / (ym(j  )-ym(j-1))

        var(39) =  areav_n(i) * (um_tau(i+1,j  )-um_tau(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um_tau(i+1,j-1)-um_tau(i,j-1)) / (xm(i+1)-xm(i))
            !bulk artificial viscosity term from Ramshaw(1990)
        var(17) = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) - b_art * (var(39) + var(37))
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - var(12) * var(35) &
                                  + var(4) * var(27) + var(6) * var(29)  &
                                  + var(8) * var(31) + var(10) * var(33))  &
                                  - var(17) &
              + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
               -epsilon1(i,j)*  ( var(35)/(Re*Darcy_number) + &
              CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* var(35)*((var(26)**2.d0 +var(35)**2.d0)**0.5d0))*liga_poros(i,j)
              
!$acc end data
return
end subroutine upwind_V1

!-------------------------------------------------------------------------------------------------------
!#######################################################################################
!#######################################################################################
subroutine upwind_V2x(um_tau,vi,p,RV,InvFr2,T,i,j)
    use openacc
    use comum
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vi, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
    real(8) :: InvFr2
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: q_art, artMAX
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
        fn = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        fs = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        fe = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        fw = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)
        df = fe - fw + fn - fs
        Dn = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))
	aw = Dw + MAX(fw , 0.0d0)
        as = Ds + MAX(fs , 0.0d0)
        ae = De + MAX(0.0d0 , -fe)
        an = Dn + MAX(0.0d0 , -fn)
        ap = aw + ae + as + an + df
        v_W = vi(i-1,j  )
        v_E = vi(i+1,j  )
        v_S = vi(i  ,j-1)
        v_N = vi(i  ,j+1)
        v_P = vi(i  ,j  )
        u_P = um_tau(i  ,j  )     !!!!!POROUS NEW
            dvdydy =  areav_n(i) * ( v_N - v_P ) / (ym(j+1)-ym(j  )) &
                     -areav_s(i) * ( v_P - v_S ) / (ym(j  )-ym(j-1))

            dydudx =  areav_n(i) * (um(i+1,j  )-um(i,j  )) / (xm(i+1)-xm(i)) &
                     -areav_s(i) * (um(i+1,j-1)-um(i,j-1)) / (xm(i+1)-xm(i))
            !bulk artificial viscosity term from Ramshaw(1990)
            q_art = epsilon1(i,j)* ( p(i,j)-p(i,j-1) )/ (y(j)-y(j-1)) - b_art * (dydudx + dvdydy)
	RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * v_P &
                                  + aw * v_W + ae * v_E  &
                                  + as * v_S + an * v_N)  &
                                  - q_art &
              + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
               -epsilon1(i,j)*  ( v_P/(Re*Darcy_number) + &
              CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)

return
end subroutine upwind_V2x

!#########################################################################
!#########################################################################

subroutine solve_V(um,vm,vm_n,um_tau,vm_tau,vm_n_tau,p,T,InvFr2,residual_v,var,res_v,Rv,vi,areav_e,areav_w,areav_n,areav_s,ym,xm,Re,x,y,liga_poros,epsilon1,artDivV)
    use comum
    use openacc
    !acc routine
    !$acc routine(RESV1)
    !$acc routine(RESV2)
    !$acc routine(upwind_v1)
    !$acc routine(upwind_v2)
    !$acc routine(bcUV2)
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm, vm_n, res_v
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, vm_n_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
    REAL(8), DIMENSION(45) :: var
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: RV, vi
    REAL(8) :: residual_v, InvFr2
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: Re
    REAL(8)            :: liga_poros(imax,jmax)
    REAL(8)            :: epsilon1(imax,jmax)
    REAL(8)            :: areav_n(1:imax)            !area n de v
    REAL(8)            :: areav_s(1:imax)            !area s de v
    REAL(8)            :: areav_e(2:jmax)            !area e de v
    REAL(8)            :: areav_w(2:jmax)            !area w de v
    REAL(8)            :: artDivV(imax,jmax)


    RV = 0.d0
    res_v = 0.d0
    vi = 0.d0
    var = 0.d0

! $acc enter data copyin(vi, RV, res_v,um,vm, vm_n, um_tau,vm_tau,vm_n_tau,p,t, v_i, areav_e, areav_w,areav_s,&
! $acc&areav_n, epsilon1, re, ym, x, xm, y, b_art, liga_poros, artDivV, var, residual_v)

!$acc parallel loop collapse(2) private(i,j, var)
    DO i=3,imax-2
       DO j=3,jmax-1
          CALL RESV1(um_tau,vm_tau,p,T,InvFr2,RV,i,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,artDivV,b_art,liga_poros,var)
       ENDDO
    ENDDO
!$acc end parallel

    !$acc parallel loop private(i, j, var)
    DO i=2,imax-1
       CALL upwind_V1(um_tau,vm_tau,p,RV,InvFr2,T,i,2,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
       CALL upwind_V1(um_tau,vm_tau,p,RV,InvFr2,T,i,jmax,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
    ENDDO
    !$acc end parallel loop
    !$acc parallel loop private(i, j, var)
    DO j=2,jmax
       CALL upwind_V1(um_tau,vm_tau,p,RV,InvFr2,T,2,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
       CALL upwind_V1(um_tau,vm_tau,p,RV,InvFr2,T,imax-1,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
    ENDDO
    !$acc end parallel loop

    !$acc parallel loop collapse(2) private(i,j) 
    DO i=2,imax-1
        DO j=3,jmax-1
            !call solve_res_v(i,j,vm,vm_tau,RV,res_v)
	    res_v(i,j) =( (vm(i,j)-vm_tau(i,j)) +  RV(i,j)*dt) * dtau 
            vi(i,j) = ( vm_tau(i,j) + res_v(i,j) ) 
        ENDDO
    ENDDO
    !$acc end parallel loop


    call bcUV2(um_tau,vi)

!$acc parallel loop collapse(2) private(i, j, var)
DO i=3,imax-2
    DO j=3,jmax-1
        CALL RESV2(um_tau,vi,p,T,InvFr2,RV,i,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,artDivV,b_art,liga_poros,var)
    ENDDO
ENDDO
!$acc end parallel loop
	
    !$acc parallel loop private (i,var)
    DO i=2,imax-1
       CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,i,2,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
       CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,i,jmax,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
    ENDDO
    !$acc end parallel loop
    !$acc parallel loop private (j,var)
    DO j=2,jmax
       CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,2,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
       CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,imax-1,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
    ENDDO
    !$acc end parallel loop

    !$acc parallel loop collapse(2) private(i,j) 
    DO i=2,imax-1
        DO j=3,jmax-1
    	    res_v(i,j) =( (vm(i,j)-vm_tau(i,j)) +  RV(i,j)*dt) * dtau 
            vi(i,j) =( 0.75d0 * vm_tau(i,j) + 0.25d0 * &
                     ( vi(i,j) + res_v(i,j)) ) 
        ENDDO
    ENDDO
    !$acc end parallel loop

    !call bcUV3(um_tau,vi)
    call bcUV2(um_tau,vi)

    !$acc parallel loop collapse(2) private(i,j,var)
    DO i=3,imax-2
      DO j=3,jmax-1
          CALL RESV2(um_tau,vi,p,T,InvFr2,RV,i,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,artDivV,b_art,liga_poros,var)
      ENDDO
    ENDDO
   !$acc end parallel loop
   !$acc parallel loop private (i,var)
    DO i=2,imax-1
       CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,i,2,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
       CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,i,jmax,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
    ENDDO
   !$acc end parallel loop
   !$acc parallel loop private (j,var)
    DO j=2,jmax
       CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,2,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
       CALL upwind_V2(um_tau,vi,p,RV,InvFr2,T,imax-1,j,areav_n,areav_s,areav_e,areav_w,epsilon1,Re,ym,x,xm,y,b_art,liga_poros,var)
    ENDDO
   !$acc end parallel loop

    !$acc parallel loop collapse(2) private(i,j,var) 
    DO i=2,imax-1
        DO j=3,jmax-1
            res_v(i,j) =( (vm(i,j)-vm_tau(i,j)) +  RV(i,j)*dt) * dtau    
            vm_n_tau(i,j) = 1.0d0 / 3.0d0 * vm_tau(i,j) + 2.0d0 / 3.0d0 * &
                     ( vi(i,j) + res_v(i,j)) 
        ENDDO
    ENDDO
    !$acc end parallel loop
   
    call bcUV2(um_tau,vm_n_tau)

!$acc parallel
    residual_v =  MAXVAL(ABS(res_v)) 
!$acc end parallel

! $acc exit data copyout(vi, RV, res_v,um,vm, vm_n, um_tau,vm_tau,vm_n_tau,p,t, v_i, areav_e, areav_w,areav_s,&
! $acc&areav_n, epsilon1, re, ym, x, xm, y, b_art, liga_poros, artDivV, var, residual_v)

! $delete(vi, RV, res_v,um,vm, vm_n, um_tau,vm_tau,vm_n_tau,p,t, v_i, areav_e, areav_w,areav_s,&
! $delete&areav_n, epsilon1, re, ym, x, xm, y, b_art, liga_poros, artDivV, var, residual_v)

    !open (550,file='data/resu.dat')
    !     do i = 2,imax-1
    !         do j = 3,jmax-1
             !write (550,*) X(i) , Y(j) , abs(res_v(i,j))
    !         enddo
    !     enddo    
    !close(550)
RETURN
END SUBROUTINE solve_V

subroutine solve_res_v(i,j,m,m_tau,RM,res_m)
    use comum
    implicit none
    integer :: i, j
   
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: m,m_tau,RM
    REAL(8), DIMENSION(2:imax-1,3:jmax-1) :: res_m
           
    res_m(i,j) =( (m(i,j)-m_tau(i,j)) +  RM(i,j)*dt) * dtau 
  
RETURN
END SUBROUTINE solve_res_v


!#######################################################################################
!#######################################################################################
!#######################################################################################
       
function alpha(x)
    !$acc routine
    implicit none
    real(8) :: alpha, x

    if(x.GT.0.0d0) then
       alpha = 1.0d0
    elseif(x.LT.0.0d0) then
       alpha = 0.0d0
    endif

return
end function alpha

!#######################################################################################
!#######################################################################################
!#######################################################################################


SUBROUTINE quick(ap,aw,aww,ae,aee,an,ann,as,ass,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
   USE COMUM
   IMPLICIT NONE
    REAL(8) :: fw, fe, fs, fn, df
    REAL(8) :: aw, aww, ae, aee, as, ass, an, ann, ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: alpha

!!! common coefficient in 3rd-order upwind QUICK Scheme
        aw = Dw + 0.75d0  * alpha(fw) * fw &
                + 0.125d0 * alpha(fe) * fe &
                + 0.375d0 * ( 1.0d0 - alpha(fw) ) * fw

        ae = De - 0.375d0* alpha(fe) * fe &
                - 0.75d0  * ( 1.0d0 - alpha(fe) ) * fe &
                - 0.125d0 * ( 1.0d0 - alpha(fw) ) * fw                

        as = Ds + 0.75d0  * alpha(fs) * fs &
                + 0.125d0 * alpha(fn) * fn &
                + 0.375d0 * ( 1.0d0 - alpha(fs) ) * fs

        an = Dn - 0.375d0* alpha(fn) * fn &
                - 0.75d0  * ( 1.0d0 - alpha(fn) ) * fn & 
                - 0.125d0 * ( 1.0d0 - alpha(fs) ) * fs

        aww = -0.125d0 *           alpha(fw)   * fw
        aee =  0.125d0 * ( 1.0d0 - alpha(fe) ) * fe
        ass = -0.125d0 *           alpha(fs)   * fs
        ann =  0.125d0 * ( 1.0d0 - alpha(fn) ) * fn

        ap = aw + ae + as + an + aww + aee + ass + ann + df

RETURN
END SUBROUTINE quick

SUBROUTINE upwind(ap,aw,ae,an,as,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
   USE COMUM
   IMPLICIT NONE
    REAL(8) :: fw, fe, fs, fn, df
    REAL(8) :: aw, ae, as, an, ap
    REAL(8) :: Dn, Ds, De, Dw

        aw = Dw + MAX(fw , 0.0d0)
        as = Ds + MAX(fs , 0.0d0)
    
        ae = De + MAX(0.0d0 , -fe)
        an = Dn + MAX(0.0d0 , -fn) 
    
    
        ap = aw + ae + as + an + df

RETURN
END SUBROUTINE upwind

!#######################################################################################
!#######################################################################################
!#######################################################################################

!!! Solve Continuity Equation
SUBROUTINE solve_P(c2,p,um_n,vm_n,pn,residual_p,RP,Pi,res_p,areau_e, areau_w, areav_n, areav_s,max_vel)
    USE comum
    USE openacc
    !acc routine
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n
    REAL(8), DIMENSION(1:imax,1:jmax) :: P, Pn
    !REAL(8) :: dudx, dvdy
    REAL(8), DIMENSION(1:imax,1:jmax)     ::  RP, Pi, res_p
    !REAL(8), DIMENSION(1:imax,1:jmax) :: k1, k2, k3, k4
    REAL(8) :: c, c2, residual_p, max_vel
    REAL(8)            :: areau_e(1:jmax)            !area e de u
    REAL(8)            :: areau_w(1:jmax)            !area w de u
    REAL(8)            :: areav_n(1:imax)            !area n de v
    REAL(8)            :: areav_s(1:imax)            !area s de v
    
    
! $acc enter data copyin(um_n, vm_n, varP, RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p)
! $acc data present(um_n, vm_n, varP,RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p)

    RP = 0.d0
    Pi = 0.d0
    res_p = 0.d0
    max_vel = MAXVAL(um_n**2.d0) + MAXVAL(vm_n**2.d0)
    ! from AN ARTIFICIAL COMPRESSIBILITY METHOD FOR INCOMPRESSIBLE FLOWS - M. M. Rahman, T. Siikonen
    c2 = beta

!!! RALSTON'S METHOD (Second Order Runge-Kutta)

    !CALL RESP(um_n,vm_n,p,RP)
        !RESP(um_n,vm_n,p,RP)
    !varP(1) = dudx
    !varP(2) = dvdy
    !$acc parallel loop collapse(2)  
    do i=2,imax-1
        do j=2,jmax-1
            !varP(1) = um_n(i+1,j) * areau_e(j) - um_n(i,j) * areau_w(j)
            !varP(2) = vm_n(i,j+1) * areav_n(i) - vm_n(i,j) * areav_s(i)
            RP(i,j) = - ( (um_n(i+1,j) * areau_e(j) - um_n(i,j) * areau_w(j)) + (vm_n(i,j+1) * areav_n(i) - vm_n(i,j) * areav_s(i)) )
         enddo
     enddo
     !$acc end parallel

    !$acc parallel loop collapse(2) 
    DO i=2,imax-1
        DO j=2,jmax-1
            pi(i,j) = p(i,j) + dtau * RP(i,j)* c2 
        ENDDO
    ENDDO
    !$acc end parallel loop

    !CALL bcP1(pi)
!--------contorno inferior e superior --------
    !$acc parallel loop private(i,j) 
    do i=1,imax
        pi(i,1  ) = pi(i,2)
        pi(i,jmax) = pi(i,jmax-1) + 1.d0*(pi(i,jmax-1)-pi(i,jmax-2))
    enddo
    !$acc end parallel loop
!-------contorno esquerdo e direito -------
   !$acc parallel loop private(i,j) 
   do j=1,jmax
        pi(1  ,j) =  pi(2,j)
        pi(imax,j) = pi(imax-1,j)
   enddo
   !$acc end parallel loop

    !CALL RESP(um_n,vm_n,pi,RP)
          !RESP(um_n,vm_n,p,RP)
    !$acc parallel loop collapse(2) private(i,j) 
    do i=2,imax-1
        do j=2,jmax-1
            !dudx = um_n(i+1,j) * areau_e(j) - um_n(i,j) * areau_w(j)
            !dvdy = vm_n(i,j+1) * areav_n(i) - vm_n(i,j) * areav_s(i)
            RP(i,j) = - ( (um_n(i+1,j) * areau_e(j) - um_n(i,j) * areau_w(j)) + (vm_n(i,j+1) * areav_n(i) - vm_n(i,j) * areav_s(i)))
         enddo
     enddo
     !$acc end parallel loop

    !$acc parallel loop collapse(2)
    DO i=2,imax-1
        DO j=2,jmax-1
            pi(i,j) = 0.75d0 * p(i,j) + 0.25d0 * (pi(i,j) + dtau * RP(i,j)* c2)
        ENDDO
    ENDDO
    !$acc end parallel loop

    !CALL bcP1(pi)
!--------contorno inferior e superior --------
   !$acc parallel loop
   do i=1,imax
        pi(i,1  ) = pi(i,2)
        pi(i,jmax) = pi(i,jmax-1) + 1.d0*(pi(i,jmax-1)-pi(i,jmax-2))
   enddo
   !$acc end parallel
!-------contorno esquerdo e direito -------
   
   !$acc parallel loop 
   do j=1,jmax
        pi(1  ,j) =  pi(2,j)
        pi(imax,j) = pi(imax-1,j)
   enddo
   !$acc end parallel
    
!CALL RESP(um_n,vm_n,pi,RP)
    !$acc parallel loop collapse(2) 
    do i=2,imax-1
        do j=2,jmax-1
            !dudx = um_n(i+1,j) * areau_e(j) - um_n(i,j) * areau_w(j)
            !dvdy = vm_n(i,j+1) * areav_n(i) - vm_n(i,j) * areav_s(i)
            RP(i,j) = - ( (um_n(i+1,j) * areau_e(j) - um_n(i,j) * areau_w(j)) + (vm_n(i,j+1) * areav_n(i) - vm_n(i,j) * areav_s(i)) )
         enddo
     enddo
     !$acc end parallel loop

    !$acc parallel loop 
    DO i=2,imax-1
        DO j=2,jmax-1
            res_p(i,j) = dtau * RP(i,j) * c2
            pn(i,j) = 1.0d0 / 3.0d0 * p(i,j) + 2.0d0 / 3.0d0 * (pi(i,j) + res_p(i,j))
        ENDDO
    ENDDO
    !$acc end parallel loop

    !CALL bcP(pn)
!--------contorno inferior e superior --------
   !$acc parallel loop 
   do i=1,imax
        pn(i,1  ) = pn(i,2)
        pn(i,jmax) = pn(i,jmax-1) + 1.d0*(pn(i,jmax-1)-pn(i,jmax-2))
   enddo
   !$acc end parallel
!-------contorno esquerdo e direito -------
   !$acc parallel loop  
   do j=1,jmax
        pn(1  ,j) =  pn(2,j)
        pn(imax,j) = pn(imax-1,j)
    enddo
   !$acc end parallel

    !$acc parallel
    residual_p = MAXVAL(ABS(res_p))
    !$acc end parallel

! $acc end data
! $acc exit data copyout(um_n, vm_n, RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p)
! $delete(um_n, vm_n, RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p)
RETURN
END SUBROUTINE solve_P


!#######################################################################################
!#######################################################################################
!#######################################################################################


SUBROUTINE RESP(um_n,vm_n,p,RP)
    USE comum
    USE openacc
    !USE omp_lib
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n
    REAL(8), DIMENSION(1:imax,1:jmax) :: P, RP !, Pi
    REAL(8) :: dudx, dvdy

    !$omp parallel do private(i,j,dudx,dvdy) 
    do i=2,imax-1
        do j=2,jmax-1
 
            dudx = um_n(i+1,j) * areau_e(j) - um_n(i,j) * areau_w(j)

            dvdy = vm_n(i,j+1) * areav_n(i) - vm_n(i,j) * areav_s(i)

            RP(i,j) = - ( dudx + dvdy )
            
         enddo
     enddo
     !$omp end parallel do

RETURN
END SUBROUTINE RESP

!#######################################################################################
!#######################################################################################
!#######################################################################################

!!! Solve Mixture Fraction
SUBROUTINE solve_Z(um_n,vm_n,Z,Z_n_tau,Z_tau, varZ,res_Z,RZ,zi,areau_e,areau_w,areav_n,areav_s,ym,xm,Pe,x,y,liga_poros,dtau)
    USE comum
    use openacc
    !acc routine
    !$acc routine(bcZ2)
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n
    REAL(8), DIMENSION(1:imax,1:jmax) :: Z!, RZ, Zi
    REAL(8), DIMENSION(1:imax,1:jmax) :: Z_n_tau, Z_tau
    !REAL(8), DIMENSION(2:imax-1,2:jmax-1) :: res_Z
    !REAL(8) :: dZudx, dZvdy
    !REAL(8) :: dZdx2, dZdy2
    !REAL(8) :: Dw,De,Ds,Dn, Dp
    !REAL(8) :: Zw,Ze,Zs,Zn, Zp
    REAL (8) :: varZ(15)
    REAL(8), DIMENSION(2:imax-1,2:jmax-1) :: res_Z
    REAL(8), DIMENSION(1:imax,1:jmax)     ::  RZ, zi
    REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    REAL(8)            :: areau_e(1:jmax)            !area e de u
    REAL(8)            :: areau_w(1:jmax)            !area w de u
    REAL(8)            :: areav_n(1:imax)            !area n de v
    REAL(8)            :: areav_s(1:imax)            !area s de v
    REAL(8)            :: Pe,dtau
    REAL(8)            :: liga_poros(imax,jmax)

 !   REAL(8), DIMENSION(1:imax,1:jmax ) :: epsilon1
!!! RALSTON'S METHOD (Second Order Runge-Kutta)
    RZ = 0.d0
    res_z = 0.d0
    zi = 0.d0

! $acc enter data copyin(dZudx, dZvdy, Z_tau, um_n, areau_e, areau_w, areau_w, areav_n, areav_s, ym, xm, Pe, x, y,Dw,De,Ds,Dn, Dp,&
! $acc&Zw,Ze,Zs,Zn, Zp, RZ, zi, T_tau, liga_poros, res_z, T, dtau, um_n_tau, vm_n_tau, T_n_tau, dx, dy)

    !CALL RESZ(um_n,vm_n,Z_tau,RZ)
!$acc parallel loop collapse(2) private(i,j,varZ) 
    do i=2,imax-1
        do j=2,jmax-1
            varZ(1) = 0.5d0 * (Z_tau(i+1,j)+Z_tau(i,j)) * um_n(i+1,j) * areau_e(j) &
                  - 0.5d0 * (Z_tau(i-1,j)+Z_tau(i,j)) * um_n(i  ,j) * areau_w(j)
            varZ(2) = 0.5d0 * (Z_tau(i,j+1)+Z_tau(i,j)) * vm_n(i,j+1) * areav_n(i) &
                  - 0.5d0 * (Z_tau(i,j-1)+Z_tau(i,j)) * vm_n(i,j  ) * areav_s(i)
            varZ(4) = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i+1)-x(i  ))
            varZ(3) = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i  )-x(i-1))
            varZ(6) = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j+1)-y(j  ))
            varZ(5) = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j  )-y(j-1))
            varZ(9) = Z_tau(i+1,j)
            varZ(8) = Z_tau(i-1,j)
            varZ(11) = Z_tau(i,j+1)
            varZ(10) = Z_tau(i,j-1)
            varZ(12) = Z_tau(i,j  )
            varZ(7) = varZ(4) + varZ(3) + varZ(6) + varZ(5)
            RZ(i,j) = 1.d0 / (xm(i+1)-xm(i)) / (ym(j+1)-ym(j)) *&
                    (-varZ(7)*varZ(12) + varZ(4)*varZ(9) + varZ(3)*varZ(8) + varZ(6)*varZ(11) + varZ(5)*varZ(10) -&
                     (1.d0-liga_poros(i,j))* (varZ(1) + varZ(2)) )!/epsilon1(i,j))  

         enddo
     enddo
     !$acc end parallel loop

    !$acc parallel loop collapse(2) private(i,j) 
    DO i=2,imax-1
        DO j=2,jmax-1
            !call solve_res_Z(i,j,Z,Z_tau,RZ,res_Z)
	    res_Z(i,j) =( (Z(i,j)- Z_tau(i,j) ) +  RZ(i,j)*dt) * dtau
            Zi(i,j) = Z_tau(i,j) + res_Z(i,j)
        ENDDO
    ENDDO
    !$acc end parallel loop

    CALL bcZ2(Zi)

    !CALL RESZ(um_n,vm_n,Zi,RZ)
    !$acc parallel loop collapse(2) private(i,j,varZ) 
    do i=2,imax-1
        do j=2,jmax-1
            varZ(1) = 0.5d0 * (Zi(i+1,j)+Zi(i,j)) * um_n(i+1,j) * areau_e(j) &
                  - 0.5d0 * (Zi(i-1,j)+Zi(i,j)) * um_n(i  ,j) * areau_w(j)
            varZ(2) = 0.5d0 * (Zi(i,j+1)+Zi(i,j)) * vm_n(i,j+1) * areav_n(i) &
                  - 0.5d0 * (Zi(i,j-1)+Zi(i,j)) * vm_n(i,j  ) * areav_s(i)
            varZ(4) = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i+1)-x(i  ))
            varZ(3) = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i  )-x(i-1))
            varZ(6) = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j+1)-y(j  ))
            varZ(5) = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j  )-y(j-1))
            varZ(9) = Zi(i+1,j)
            varZ(8) = Zi(i-1,j)
            varZ(11) = Zi(i,j+1)
            varZ(10) = Zi(i,j-1)
            varZ(12) = Zi(i,j  )
            varZ(7) = varZ(4) + varZ(3) + varZ(6) + varZ(5)
            RZ(i,j) = 1.d0 / (xm(i+1)-xm(i)) / (ym(j+1)-ym(j)) *&
                    (-varZ(7)*varZ(12) + varZ(4)*varZ(9) + varZ(3)*varZ(8) + varZ(6)*varZ(11) + varZ(5)*varZ(10) -&
                     (1.d0-liga_poros(i,j))* (varZ(1) + varZ(2)) )!/epsilon1(i,j))    

         enddo
     enddo
     !$acc end parallel loop

    !$acc parallel loop collapse(2) private(i,j) 
    DO i=2,imax-1
        DO j=2,jmax-1       
            res_Z(i,j) =( (Z(i,j)-Z_tau(i,j)) +  RZ(i,j)*dt) * dtau
            Zi(i,j) = 0.75d0 * Z_tau(i,j) + 0.25d0 * (Zi(i,j) + res_Z(i,j))
        ENDDO
    ENDDO
    !$acc end parallel loop

    CALL bcZ2(Zi)

    !CALL RESZ(um_n,vm_n,Zi,RZ)
    !$acc parallel loop collapse(2) private(i,j,varZ) 
    do i=2,imax-1
        do j=2,jmax-1
            varZ(1) = 0.5d0 *(Zi(i+1,j)+Zi(i,j)) * um_n(i+1,j) * areau_e(j) - 0.5d0  *(Zi(i-1,j)+Zi(i,j)) * um_n(i  ,j) * areau_w(j)
            varZ(2) = 0.5d0 *(Zi(i,j+1)+Zi(i,j)) * vm_n(i,j+1) * areav_n(i) - 0.5d0 *(Zi(i,j-1)+Zi(i,j)) * vm_n(i,j  ) * areav_s(i)
            varZ(4) = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i+1)-x(i  ))
            varZ(3) = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i  )-x(i-1))
            varZ(6) = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j+1)-y(j  ))
            varZ(5) = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j  )-y(j-1))
            varZ(9) = Zi(i+1,j)
            varZ(8) = Zi(i-1,j)
            varZ(11) = Zi(i,j+1)
            varZ(10) = Zi(i,j-1)
            varZ(12) = Zi(i,j  )
            varZ(7) = varZ(4) + varZ(3) + varZ(6) + varZ(5)
            RZ(i,j) = 1.d0 / (xm(i+1)-xm(i)) / (ym(j+1)-ym(j)) *&
                    (-varZ(7)*varZ(12) + varZ(4)*varZ(9) + varZ(3)*varZ(8) + varZ(6)*varZ(11) + varZ(5)*varZ(10) -&
                     (1.d0-liga_poros(i,j))* (varZ(1) + varZ(2)) )!/epsilon1(i,j))  

         enddo
     enddo
     !$acc end parallel loop

    !$acc parallel loop private(i,j) 
    DO i=2,imax-1
        DO j=2,jmax-1            
            res_Z(i,j) =( (Z(i,j)-Z_tau(i,j)) +  RZ(i,j)*dt) * dtau
            Z_n_tau(i,j) = 1.0d0 / 3.0d0 * Z_tau(i,j) + 2.0d0 / 3.0d0 * (Zi(i,j) + res_Z(i,j))
        ENDDO
    ENDDO
    !$acc end parallel loop

    CALL bcZ1(Z_n_tau)

!  $ acc exit data copyout(dZudx, dZvdy, Z_tau, um_n, areau_e, areau_w, areau_w, areav_n, areav_s, ym, xm, Pe, x, y,Dw,De,Ds,Dn, Dp,&
!  $ acc&Zw,Ze,Zs,Zn, Zp, RZ, zi, T_tau, liga_poros, res_z, T, dtau, um_n_tau, vm_n_tau, T_n_tau, dx, dy)


!  $ delete(dZudx, dZvdy, Z_tau, um_n, areau_e, areau_w, areau_w, areav_n, areav_s, ym, xm, Pe, x, y,Dw,De,Ds,Dn, Dp,&
!  $  delete&Zw,Ze,Zs,Zn, Zp, RZ, zi, T_tau, liga_poros, res_z, T, dtau, um_n_tau, vm_n_tau, T_n_tau, dx, dy)

RETURN
END SUBROUTINE solve_Z

subroutine solve_res_Z(i,j,m,m_tau,RM,res_m)


    use comum
    implicit none
    integer :: i, j
   
    REAL(8), DIMENSION(1:imax  ,1:jmax) :: m,m_tau,RM
    REAL(8), DIMENSION(2:imax-1,2:jmax-1) :: res_m
            
    !res_m(i,j) =( (m(i,j)-m_tau(i,j))/dt +  RM(i,j)) * dtau 
    res_m(i,j) =( (m(i,j)-m_tau(i,j)) +  RM(i,j)*dt) * dtau 
    !res_m(i,j) = RM(i,j) * dt 

RETURN
END SUBROUTINE solve_res_Z

!#######################################################################################
!#######################################################################################
!#######################################################################################

SUBROUTINE RESZ(um_n,vm_n,Z,RZ)
    USE comum
    USE openacc
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n
    REAL(8), DIMENSION(1:imax,1:jmax) :: Z, RZ !, Zi
 !   REAL(8), DIMENSION(0:imax+1,0:jmax+1) ::epsilon1
    REAL(8) :: dZudx, dZvdy
    REAL(8) :: dZdx2, dZdy2
    REAL(8) :: Dw,De,Ds,Dn, Dp
    REAL(8) :: Zw,Ze,Zs,Zn, Zp

    !$omp parallel do private(i,j, dZudx, dZvdy, De, Dw, Dn, Ds, Ze,Zw,Zn,Zs,Zp, Dp) 
    do i=2,imax-1
        do j=2,jmax-1
            dZudx = 0.5d0 * (Z(i+1,j)+Z(i,j)) * um_n(i+1,j) * areau_e(j) &
                  - 0.5d0 * (Z(i-1,j)+Z(i,j)) * um_n(i  ,j) * areau_w(j)
            dZvdy = 0.5d0 * (Z(i,j+1)+Z(i,j)) * vm_n(i,j+1) * areav_n(i) &
                  - 0.5d0 * (Z(i,j-1)+Z(i,j)) * vm_n(i,j  ) * areav_s(i)
            De = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i+1)-x(i  ))  
            Dw = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i  )-x(i-1))  
            Dn = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j+1)-y(j  ))  
            Ds = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j  )-y(j-1))  
            Ze = Z(i+1,j)
            Zw = Z(i-1,j)
            Zn = Z(i,j+1)
            Zs = Z(i,j-1) 
            Zp = Z(i,j  ) 
            Dp = De + Dw + Dn + Ds 
            RZ(i,j) = 1.d0 / (xm(i+1)-xm(i)) / (ym(j+1)-ym(j)) *&
                    (-Dp*Zp + De*Ze + Dw*Zw + Dn*Zn + Ds*Zs - (1.d0-liga_poros(i,j))* (dZudx + dZvdy) )!/epsilon1(i,j))  

         enddo
     enddo
     !$omp end parallel do

RETURN
END SUBROUTINE RESZ


!#######################################################################################
!#######################################################################################
!#######################################################################################

SUBROUTINE solve_H(H,Z)
    USE comum
    use openacc
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: Z, H
    REAL(8) :: H_o, H_f

    !calculo de H utilizando Z
    !só vale para caso com Le = 1

    H_o = (1.d0 / q ) * Tinf + 1.d0 / (S+1.d0)

    H_f = ( 1.d0 / q ) * Tsup + 1.d0 / (S+1.d0)

    H = H_o + ( (H_f - H_o) * Z )

    call bcH(H)

RETURN
END SUBROUTINE solve_H


!#######################################################################################
!#######################################################################################
!#######################################################################################


SUBROUTINE comp_T(Z,H,T)
    USE comum
    use openacc
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax,1:jmax) :: Z, H, T

    !#omp parallel do private(i,j)
    do i=1,imax
        do j=1,jmax
            if (Z(i,j).gt.(1.d0/(S+1))) then

            T(i,j) = q * ( H(i,j) - (Z(i,j)/S) + 1.d0 / (S*(S+1.d0)) )

            else

            T(i,j) = q * ( H(i,j) + Z(i,j) - 1.d0 / (S+1.d0) )

            endif

        enddo
    enddo
    !#omp end parallel do

RETURN
END SUBROUTINE comp_T
