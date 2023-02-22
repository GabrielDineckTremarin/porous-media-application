subroutine RESU(um,vm,p,RU)

    use comum
    use omp_lib
    implicit none
    integer :: i, j
   
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: P
!    REAL(8), DIMENSION(1:imax+1,1:jmax+1) ::epsilon1,Darcy
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: q_art, artDivU(imax,jmax), artDivV(imax,jmax), artMAX
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

        CALL quick(ap,aw,aww,ae,aee,an,ann,as,ass,fw,fe,fn,fs,df,Dw,De,Dn,Ds)

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


!TIRAR
!DO i=3,imax-1
!    DO j=3,jmax-2
!    IF (flag(i,j) .EQ. C_BS) THEN
!
!        CALL upwind_U(um,vm,p,RU,i,j)
!
!    ENDIF
!    ENDDO
!ENDDO    

!mudar ordem no contorno.

!$omp parallel do private(i) 
DO i=2,imax
    CALL upwind_U(um,vm,p,RU,i,2)
    CALL upwind_U(um,vm,p,RU,i,jmax-1)
ENDDO
!$omp end parallel do

!$omp parallel do private(j) 
DO j=2,jmax-1
    CALL upwind_U(um,vm,p,RU,2,j)
    CALL upwind_U(um,vm,p,RU,imax,j)
ENDDO
!$omp end parallel do
RETURN
END SUBROUTINE RESU

!#######################################################################################
!#######################################################################################
!#######################################################################################

subroutine upwind_U(um,vm,p,RU,i,j)

    use comum
    implicit none
    integer :: i, j, itc
    
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax,1:jmax) :: P
 !   REAL(8), DIMENSION(1:imax+1,1:jmax+1) ::epsilon1,Darcy
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: alpha
    real(8) :: q_art, artDivU(imax,jmax), artDivV(imax,jmax), artMAX
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


        CALL upwind(ap,aw,ae,an,as,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
    
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
!  

  
        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * u_P &
                                  + aw * u_W + ae * u_E  &
                                  + as * u_S + an * u_N) &
                                  - q_art &
               - epsilon1(i,j)* (  u_P/(Re*Darcy_number) +&
        CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* u_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
RETURN
END SUBROUTINE upwind_U


!#######################################################################################
!#######################################################################################
!#######################################################################################

subroutine solve_U(um,vm,um_n,um_tau,vm_tau,um_n_tau,p,residual_u)


    use comum
    use omp_lib
    implicit none
    integer :: i, j
   
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um,um_n,ui, RU
    REAL(8), DIMENSION(3:imax-1,2:jmax-1) :: res_u
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax  ) :: p
!    REAL(8), DIMENSION(1:imax+1,1:jmax+1) ::epsilon1,Darcy
    REAL(8) :: residual_u

    RU = 0.d0
    res_u = 0.d0



! EXCLUIR E CRIAR UMA MATRIZ porosidade(i,j) Darcy_term(i,j)
!------ exclude body from residuals -----------
!	DO i=3,imax-1
! 	   DO j=2,jmax-1
! 	   IF (flag(i,j) .NE. C_F) THEN
!    
! 	      epsilon1= porosidade
! 	      Darcy_term = 1.d0
! 	      ELSE
! 	      epsilon1 = 1.d0
! 	       Darcy_term = 0.d0
! 	   ENDIF!
!  	  ENDDO
!       ENDDO	
!!!!!!!!!!!!!!!!!!!!!!!!!!1excluir       
       
         
    CALL RESU(um_tau,vm_tau,p,RU)


    !$omp parallel do private(i,j) 
    DO i=3,imax-1
        DO j=2,jmax-1
            
            call solve_res_u(i,j,um,um_tau,RU,res_u)
            
            ui(i,j) =  ( um_tau(i,j) + res_u(i,j))
 !IF (flag(i,j) .NE. C_F) THEN   ! AQUI BODY
 !	     ui(i,j) = 0.D0   ! AQUI BODY
!
!		ENDIF  


        ENDDO
    ENDDO
    !$omp end parallel do

    call bcUV(ui,vm_tau)

    CALL RESU(ui,vm_tau,p,RU)

    !$omp parallel do private(i,j) 
    DO i=3,imax-1
        DO j=2,jmax-1
    
            call solve_res_u(i,j,um,um_tau,RU,res_u)
                
            ui(i,j) = (0.75d0 * um_tau(i,j) + 0.25d0 * &
                      ( ui(i,j) + res_u(i,j) ) )
! IF (flag(i,j) .NE. C_F) THEN   ! AQUI BODY
! 	     ui(i,j) = 0.D0   ! AQUI BODY
!
!		ENDIF  

        ENDDO
    ENDDO
    !$omp end parallel do


    call bcUV(ui,vm_tau)

    CALL RESU(ui,vm_tau,p,RU)

    !$omp parallel do private(i,j) 
    DO i=3,imax-1
        DO j=2,jmax-1

            call solve_res_u(i,j,um,um_tau,RU,res_u)

            um_n_tau(i,j) = (1.0d0 / 3.0d0 * um_tau(i,j) + 2.0d0 / 3.0d0 * &
                      ( ui(i,j) + res_u(i,j) ) )
!      IF (flag(i,j) .NE. C_F) THEN   ! AQUI BODY
! 	     um_n_tau(i,j) = 0.D0   ! AQUI BODY
!
!		ENDIF  

  			 ENDDO
    ENDDO
    !$omp end parallel do

    call bcUV(um_n_tau,vm_tau)
   

    residual_u =  MAXVAL(ABS(res_u)) 



!    open (550,file='data/resu.dat')
!        
!         do i = 1,imax
!             do j = 1,jmax
!            
!             write (550,*) X(i) , Y(j) , res_u(i,j)
!
!             enddo
!         enddo    
!
!    close(550)

RETURN
END SUBROUTINE solve_U

subroutine solve_res_u(i,j,m,m_tau,RM,res_m)


    use comum
    implicit none
    integer :: i, j
   
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: m,m_tau,RM
    REAL(8), DIMENSION(3:imax-1,2:jmax-1) :: res_m
            
    !res_m(i,j) =( (m(i,j)-m_tau(i,j))/dt +  RM(i,j)) * dtau 
    res_m(i,j) =( (m(i,j)-m_tau(i,j)) +  RM(i,j)*dt) * dtau 
    !res_m(i,j) =( RM(i,j)) * dtau 
!    IF (flag(i,j) .NE. C_F) THEN   ! AQUI BODY
! 	     res_m(i,j) = 0.D0   ! AQUI BODY
!
!		ENDIF  

RETURN
END SUBROUTINE solve_res_u
!#######################################################################################
!#######################################################################################
!#######################################################################################



subroutine RESV(um,vm,p,T,InvFr2,RV)


    use comum
    use omp_lib
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
    real(8) :: q_art, artDivU(imax,jmax), artDivV(imax,jmax), artMAX
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx


!$omp parallel do private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dvdydy,dydudx,q_art,&
!$omp&v_W,v_WW,v_E,v_EE,v_N,v_NN,v_S,v_SS,v_P,u_P,ap,aw,aww,ae,aee,an,ann,as,ass)
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

        CALL quick(ap,aw,aww,ae,aee,an,ann,as,ass,fw,fe,fn,fs,df,Dw,De,Dn,Ds)

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

!@EXCLUIR !!!!!!!!!!!!
!DO i=3,imax-2
!    DO j=3,jmax-1
!
!    IF (flag(i,j) .EQ. C_BS) THEN!
!
!        CALL upwind_V(um,vm,p,RV,InvFr2,T,i,j)
!
!    ENDIF
!    ENDDO
!ENDDO    

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
!#######################################################################################
!#######################################################################################


subroutine upwind_V(um,vm,p,RV,InvFr2,T,i,j)

    use comum
    implicit none
    integer :: i, j, itc
    
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
    real(8) :: q_art, artDivU(imax,jmax), artDivV(imax,jmax), artMAX
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


        CALL upwind(ap,aw,ae,an,as,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
 
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
!#######################################################################################
!#######################################################################################

subroutine solve_V(um,vm,vm_n,um_tau,vm_tau,vm_n_tau,p,T,InvFr2,residual_v)


    use comum
    use omp_lib
    implicit none
    integer :: i, j
   
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm, vm_n, vi, RV
    REAL(8), DIMENSION(2:imax-1,3:jmax-1) :: res_v
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, vm_n_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: P,T
!    REAL(8), DIMENSION(1:imax+1  ,1:jmax+1) :: epsilon1,Darcy
    REAL(8) :: residual_v
    REAL(8) :: InvFr2

    RV = 0.d0
    res_v = 0.d0

!------ exclude body from residuals -----------
!EXCLUIR cria matriz

!DO i=2,imax-1
!    DO j=3,jmax-1
!    IF (flag(i,j) .NE. C_F) THEN
!
!        !res_v(i,j)=0.d0
!        epsilon1= porosidade
!             Darcy_term = 1.d0
!         ELSE
! 	      epsilon1 = 1.d0
!     Darcy_term = 0.d0
!
!    ENDIF
!    ENDDO
!ENDDO    
    CALL RESV(um_tau,vm_tau,p,T,InvFr2,RV)


    !$omp parallel do private(i,j) 
    DO i=2,imax-1
        DO j=3,jmax-1
            call solve_res_v(i,j,vm,vm_tau,RV,res_v)

            vi(i,j) = ( vm_tau(i,j) + res_v(i,j) ) 
!IF (flag(i,j) .NE. C_F) THEN   ! AQUI BODY
 !	     vi(i,j) = 0.D0   ! AQUI BODY
!
!		ENDIF  

        ENDDO
    ENDDO
    !$omp end parallel do

    call bcUV(um_tau,vi)

    CALL RESV(um_tau,vi,p,T,InvFr2,RV)

    !$omp parallel do private(i,j) 
    DO i=2,imax-1
        DO j=3,jmax-1
            call solve_res_v(i,j,vm,vm_tau,RV,res_v)
    
            vi(i,j) =( 0.75d0 * vm_tau(i,j) + 0.25d0 * &
                     ( vi(i,j) + res_v(i,j)) )
!     IF (flag(i,j) .NE. C_F) THEN   ! AQUI BODY
! 	     vi(i,j) = 0.D0   ! AQUI BODY
!
!		ENDIF  

        ENDDO
    ENDDO
    !$omp end parallel do


    call bcUV(um_tau,vi)

    CALL RESV(um_tau,vi,p,T,InvFr2,RV)

    !$omp parallel do private(i,j) 
    DO i=2,imax-1
        DO j=3,jmax-1
            call solve_res_v(i,j,vm,vm_tau,RV,res_v)
                
            vm_n_tau(i,j) = 1.0d0 / 3.0d0 * vm_tau(i,j) + 2.0d0 / 3.0d0 * &
                     ( vi(i,j) + res_v(i,j)) 
!     IF (flag(i,j) .NE. C_F) THEN   ! AQUI BODY
! 	     vm_n_tau(i,j) = 0.D0   ! AQUI BODY
!
!		ENDIF  

        ENDDO
    ENDDO
    !$omp end parallel do
   
call bcUV(um_tau,vm_n_tau)



residual_v =  MAXVAL(ABS(res_v)) 

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
            
    !res_m(i,j) =( (m(i,j)-m_tau(i,j))/dt +  RM(i,j)) * dtau 
    res_m(i,j) =( (m(i,j)-m_tau(i,j)) +  RM(i,j)*dt) * dtau 
    !res_m(i,j) =( RM(i,j)) * dtau 
  !    IF (flag(i,j) .NE. C_F) THEN   ! AQUI BODY
 !	     res_m(i,j) = 0.D0   ! AQUI BODY
!
!		ENDIF  

RETURN
END SUBROUTINE solve_res_v


!#######################################################################################
!#######################################################################################
!#######################################################################################
       
function alpha(x)
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
SUBROUTINE solve_P(c2,p,um_n,vm_n,pn,residual_p)
    USE comum
    use omp_lib
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n
    REAL(8), DIMENSION(1:imax,1:jmax) :: P, Pn, RP, Pi, res_p

    REAL(8), DIMENSION(1:imax,1:jmax) :: k1, k2, k3, k4
    REAL(8) :: c, c2, residual_p, max_vel
res_p = 0.d0

    max_vel = MAXVAL(um_n**2.d0) + MAXVAL(vm_n**2.d0)

    ! from AN ARTIFICIAL COMPRESSIBILITY METHOD FOR INCOMPRESSIBLE FLOWS - M. M. Rahman, T. Siikonen
    !c = beta !* SQRT(MAX(max_vel,(0.5d0 * v_i**2.d0))) 

    !c2 = c*c
    c2 = beta !* MAXVAL(vm_n**2.d0)

!------ exclude body from residuals -----------


!DO i=2,imax-1
!    DO j=3,jmax-1
!    IF (flag(i,j) .NE. C_F) THEN
!
!       ! res_p(i,j)=0.d0
!		epsilon1= porosidade
!		     Darcy_term = 1.d0
!         ELSE
! 	      epsilon1 = 1.d0
! 	           Darcy_term = 0.d0
!    ENDIF
!    ENDDO
!ENDDO    

!!! RALSTON'S METHOD (Second Order Runge-Kutta)
!mass conservation

    CALL RESP(um_n,vm_n,p,RP)

    !$omp parallel do private(i,j) 
    DO i=2,imax-1
        DO j=2,jmax-1
            
            pi(i,j) = p(i,j) + dtau * RP(i,j)* c2 

        ENDDO
    ENDDO
    !$omp end parallel do

    CALL bcP(pi)

    CALL RESP(um_n,vm_n,pi,RP)

    !$omp parallel do private(i,j) 
    DO i=2,imax-1
        DO j=2,jmax-1
                
            pi(i,j) = 0.75d0 * p(i,j) + 0.25d0 * (pi(i,j) + dtau * RP(i,j)* c2)

        ENDDO
    ENDDO
    !$omp end parallel do

    CALL bcP(pi)

    CALL RESP(um_n,vm_n,pi,RP)
    
    !$omp parallel do private(i,j) 
    DO i=2,imax-1
        DO j=2,jmax-1
            res_p(i,j) = dtau * RP(i,j) * c2
            
            pn(i,j) = 1.0d0 / 3.0d0 * p(i,j) + 2.0d0 / 3.0d0 * (pi(i,j) + res_p(i,j))
        ENDDO
    ENDDO
    !$omp end parallel do

    CALL bcP(pn)

    residual_p = MAXVAL(ABS(res_p))

    


RETURN
END SUBROUTINE solve_P


!#######################################################################################
!#######################################################################################
!#######################################################################################


SUBROUTINE RESP(um_n,vm_n,p,RP)
    USE comum
    USE omp_lib
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n
    REAL(8), DIMENSION(1:imax,1:jmax) :: P, RP, Pi
    REAL(8) :: dudx, dvdy

    !$omp parallel do private(i,j,dudx,dvdy) 
    do i=2,imax-1
        do j=2,jmax-1
 
            dudx = um_n(i+1,j) * areau_e(j) - um_n(i,j) * areau_w(j)

            dvdy = vm_n(i,j+1) * areav_n(i) - vm_n(i,j) * areav_s(i)

!AQUII  epsilon
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
SUBROUTINE solve_Z(um_n,vm_n,Z,Z_n_tau,Z_tau)
    USE comum
    use omp_lib
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n
    REAL(8), DIMENSION(1:imax,1:jmax) :: Z, RZ, Zi
    REAL(8), DIMENSION(1:imax,1:jmax) :: Z_n_tau, Z_tau
    REAL(8), DIMENSION(2:imax-1,2:jmax-1) :: res_Z

 !   REAL(8), DIMENSION(1:imax,1:jmax ) :: epsilon1
!!! RALSTON'S METHOD (Second Order Runge-Kutta)
    CALL RESZ(um_n,vm_n,Z_tau,RZ)

    !$omp parallel do private(i,j) 
    DO i=2,imax-1
        DO j=2,jmax-1
            call solve_res_Z(i,j,Z,Z_tau,RZ,res_Z)
            Zi(i,j) = Z_tau(i,j) + res_Z(i,j)
!            solve_res_Z(i,j,m,m_tau,RM,res_m)
        ENDDO
    ENDDO
    !$omp end parallel do

    CALL bcZ(Zi)

    CALL RESZ(um_n,vm_n,Zi,RZ)

    !$omp parallel do private(i,j) 
    DO i=2,imax-1
        DO j=2,jmax-1
                
            call solve_res_Z(i,j,Z,Z_tau,RZ,res_Z)

            Zi(i,j) = 0.75d0 * Z_tau(i,j) + 0.25d0 * (Zi(i,j) + res_Z(i,j))

        ENDDO
    ENDDO
    !$omp end parallel do

    CALL bcZ(Zi)

    CALL RESZ(um_n,vm_n,Zi,RZ)
    !$omp parallel do private(i,j) 
    DO i=2,imax-1
        DO j=2,jmax-1
            
            call solve_res_Z(i,j,Z,Z_tau,RZ,res_Z)

            Z_n_tau(i,j) = 1.0d0 / 3.0d0 * Z_tau(i,j) + 2.0d0 / 3.0d0 * (Zi(i,j) + res_Z(i,j))

 !	   IF (flag(i,j) .NE. C_F) THEN
 ! 	     Z_n_tau(i,j) = Temp_cylinder
 !
 !		ENDIF
        ENDDO
    ENDDO
    !$omp end parallel do

    

    CALL bcZ(Z_n_tau)


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
    USE omp_lib
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n
    REAL(8), DIMENSION(1:imax,1:jmax) :: Z, RZ, Zi
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

!            dZudx = Z(i+1,j) * um_n(i+1,j) * areau_e(j) &
!                  - Z(i-1,j) * um_n(i  ,j) * areau_w(j)
!
!            dZvdy = Z(i,j+1) * vm_n(i,j+1) * areav_n(i) &
!                  - Z(i,j-1) * vm_n(i,j  ) * areav_s(i)

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
