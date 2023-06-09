!!!    ArtIFicial Compressibility Methods
!!!    Solve Momentum Equation with QUICK Scheme
!!!
!!!              Oxidant
!!!         |---------------|
!!!         |               |
!!!         |..      |      |
!!!         |  .     |      | 
!!!  Symm   |   .    | g    | Oxidant
!!!         |   .    v      |
!!!         |  .            |
!!!         |..cylinder     |
!!!         |               |
!!!         |---------------|
!!!              Oxidant
!!!  BASED ON:
!!!Versteeg, H. K., and W. Malalasekera. 
!!!"An introduction to computationnal Fluid Dynamics, The finite volume control, éd." (1995).


PROGRAM main
USE openacc
USE comum
!$acc routine(solve_p)
!$acc routine(solve_z)
!$acc routine(solve_v)
!$acc routine(solve_u)
!$acc routine(convergence)
!$acc routine(lacointerno)
!$acc routine(terat_convergence)

IMPLICIT NONE

    INTEGER :: tr, i, j !itc 
    INTEGER*4 today(3), now(3)
   
    !REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um, um_n, res_u
    !REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm, vm_n, res_v
    
    !REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, um_n_tau
    !REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, vm_n_tau

    !REAL(8), DIMENSION(1:imax,1:jmax)     :: u, v, P, Pn, H, T,Z
    !REAL(8), DIMENSION(1:imax,1:jmax)     :: T_n_tau, T_tau 
    !REAL(8)                               :: InvFr2

    !REAL(8) :: residual_p, residual_u,residual_v
    !REAL(8) :: c2, error
    REAL(8) :: artMAX, yf, xf

    character(len=128) :: pwd
        REAL ETIME,clockTIME,TARRAY(2)
        clockTIME=ETIME(TARRAY)


        CALL idate(today)   ! today(1)=day, (2)=month, (3)=year
        CALL itime(now)     ! now(1)=hour, (2)=minute, (3)=second
        CALL get_environment_variable('PWD',pwd)       

        CALL init

!------------------------------

        Too  = Ts * nToo
        Tsup = Ts  / Too
        Tinf = Too / Too

!------------------------------

        CALL properties    

        S =  ( nu * YF_b) / YO_oo

!!----- flame width and nat conv velocity as reference ------
!        L_c = ( S + 1.d0) * ao
!
!        v_c = DSQRT( L_c * DABS(g) )
!!-----------------------------------------------------------

!----- burner as reference -----
     !   Pe = 2.d0
       v_i = 0.7d0 !v_c / v_c
        L_c = ao
        
        v_c = v_i !Pe * alpha_tot / L_c
!------------------------------

        q = (q_dim * YF_b) / (cp_tot * Too)

        Pr = nu_tot / alpha_tot

        Re = v_c * L_c / nu_tot
        
        Pe = Re * Pr

        Fr = v_c / sqrt( Dabs(g) * L_c )

        InvFr2 = 1.d0 / (Fr*Fr)
      

        5 FORMAT ( 'Date ', i2.2, '/', i2.2, '/', i4.4, &
                   '; time ',i2.2, ':', i2.2, ':', i2.2 )

            WRITE(*,*) '----------------------'
            WRITE(*,5)  today(2), today(1), today(3), now
            WRITE(*,*) 'Current working directory: ',trim(pwd)
            WRITE(*,*) '----------------------'
            WRITE(*,*) 'Tsup =',Tsup
            WRITE(*,*) 'Tinf =',Tinf
            WRITE(*,*) '----------------------'
            WRITE(*,*) 'Mesh =',imax,'x',jmax 
            WRITE(*,*) 'imax * jmax =', int(imax*jmax)
            WRITE(*,*) 'eps =',eps
            WRITE(*,*) '----------------------'
            WRITE(*,*) 'alpha =',alpha_tot  , '[m^2/s]'
            WRITE(*,*) 'cp =   ',cp_tot     , '[J/kgK]'
            WRITE(*,*) 'k =    ',k_tot      , '[W/mK]'
            WRITE(*,*) 'rho =  ',rho_tot    , '[kg/m^3]'
            WRITE(*,*) 'nu =   ',nu_tot     , '[m^2/s]'
            WRITE(*,*) 'mu =   ',nu_tot*rho_tot, '[kg/ms]'
            WRITE(*,*) '----------------------'
            WRITE(*,*) 'v_c ='  ,v_c        , '[m/s]'
            WRITE(*,*) 'L_c ='  ,L_c        , '[m]'
            WRITE(*,*) 't_c ='  ,L_c/v_c    , '[s]'
            WRITE(*,*) 'v_i ='  ,v_i
            WRITE(*,*) 'V_idim ='  ,v_i*v_c , '[m/s]'
            WRITE(*,*) '----------------------'
            WRITE(*,*) 'g = '  ,g , '[m/s^2]'
            WRITE(*,*) 'S = '  ,S
            WRITE(*,*) 'q = '  ,q
            WRITE(*,*) 'Pr ='  ,Pr
            WRITE(*,*) 'Re ='  ,Re
            WRITE(*,*) 'Pe ='  ,Pe
            WRITE(*,*) 'Fr ='  ,Fr
            WRITE(*,*) 'InvFr^2 =',InvFr2
            WRITE(*,*) '----------------------'
        

        itc = 1 !iteração inicial
        error=100.00d0
        tr = 1
        time = 0.d0
        residual_p = 0.d0

!!--------- create mesh -------------------------------------       
        CALL mesh

DO i=1,imax
 	   DO j=1,jmax
 	   IF (flag(i,j) .NE. C_F) THEN
 	     T(i,j) = Temp_cylinder
 	   ENDIF
  	  ENDDO
	ENDDO 


!!-----------------------------------------------------------

!!-------- set up initial flow field --------------------------------------------
        IF (start_mode .EQ. 0) THEN    
            CALL IC(um,vm,p,T)
        ELSEIF (start_mode .EQ. 1) THEN
            CALL restart(um,vm,p,T)
        ELSEIF (start_mode .EQ. 2) THEN
            CALL restart_dom(um,vm,p,T,I,H)
        ENDIF
            
!------pseudo time step ---------------------------------------------------------        
!        dtau = dtau_f * minval( [MINVAL(areau_s)*MINVAL(areau_s)/(4.d0*Pe),&
!                      MINVAL(areau_s)/maxval(abs(um)) ,&
!                      MINVAL(areav_w)/maxval(abs(vm)) ])
dtau = 8.d-2

!------physical time step ---------------------------------------------------------        
!        dt = 1.d4 * minval( [MINVAL(areau_s)*MINVAL(areau_s)/(4.d0*Pe),&
!                      MINVAL(areau_s)/maxval(abs(um)) ,&
!                      MINVAL(areav_w)/maxval(abs(vm)) ])
dt = 1.d-2

!!-------- physical-time calculation starts ---------------------------------------------------
  !      CALL solve_H(H,Z)


   !     CALL comp_T(Z,H,T)

            um_tau = um
            vm_tau = vm


!$acc update device( c2, dtau, beta, itc, error, time, dt, final_time)

varM(1) = itc

!$acc enter data copyin(um, um_n, vm, vm_n, RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p,&
!$acc& Z_tau, varZ, ym, xm, Pe, x, y, RZ, zi, T_tau, liga_poros, res_z, T, dtau, um_n_tau, vm_n_tau, T_n_tau, dx, dy,&
!$acc& um_tau, vm_tau, residual_u, RU, res_u, ui, epsilon1, b_art, Re, artDivU, v_i, var,&
!$acc& vi, RV, res_v,t, v_i, Re, artDivV, var, residual_v, itc, varM)


DO WHILE (time .LT. final_time)

     time = time + dt

! $acc enter data copyin(um, um_n, vm, vm_n, RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p,&
! $acc& Z_tau, varZ, ym, xm, Pe, x, y, RZ, zi, T_tau, liga_poros, res_z, T, dtau, um_n_tau, vm_n_tau, T_n_tau, dx, dy,&
! $acc& um_tau, vm_tau, residual_u, RU, res_u, ui, epsilon1, b_art, Re, artDivU, v_i, var,&
! $acc& vi, RV, res_v,t, v_i, Re, artDivV, var, residual_v, itc, varM)

    !!-------- pseudo-time calculation starts ---------------------------------------------------
        DO WHILE(itc.LT.itc_max) !itc
          
! $acc enter data copyin(um, um_n, vm, vm_n, RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p,&
! $acc& Z_tau, varZ, ym, xm, Pe, x, y, RZ, zi, T_tau, liga_poros, res_z, T, dtau, um_n_tau, vm_n_tau, T_n_tau, dx, dy,&
! $acc& um_tau, vm_tau, residual_u, RU, res_u, ui, epsilon1, b_art, Re, artDivU, v_i, var,&
! $acc& vi, RV, res_v,t, v_i, Re, artDivV, var, residual_v, itc, varM)


!-------- Solve Momentum Equation with QUICK Scheme
     CALL solve_U(um,vm,um_n,um_tau,vm_tau,um_n_tau,pn,residual_u,var,res_u,RU,UI,areau_e,areau_w,areau_n,areau_s,ym,xm,Re,x,y,liga_poros,dtau,epsilon1,artDivU)

     CALL solve_V(um,vm,vm_n,um_tau,vm_tau,vm_n_tau,pn,T,InvFr2,residual_v,var,res_v,Rv,vi,areav_e,areav_w,areav_n,areav_s,ym,xm,Re,x,y,liga_poros,epsilon1,artDivV)


     CALL solve_P(c2,p,um_n_tau,vm_n_tau,pn,residual_p,RP,Pi,res_p,areau_e, areau_w, areav_n, areav_s, max_vel) 

    !-------- Solve Energy Equation ------------------------------------------------
     CALL solve_Z(um_n_tau,vm_n_tau,T,T_n_tau,T_tau,varZ,res_Z,RZ,zi,areau_e,areau_w,areav_n,areav_s,ym,xm,Pe,x,y,liga_poros,dtau)
 
    !    CALL solve_H(H,T_n_tau)
    !    CALL comp_T(Z_n_tau,H,T)
           
    !!-------- check convergence

    CALL convergence(itc,c2,error,residual_p,residual_u,residual_v)
   
    !CALL convergencia(varM, residual_u, residual_v, residual_p, itc)

! $acc exit data copyout(um, um_n, vm, vm_n, RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p,&
! $acc& Z_tau, varZ, ym, xm, Pe, x, y, RZ, zi, T_tau, liga_poros, res_z, T, dtau, um_n_tau, vm_n_tau, T_n_tau, dx, dy,&
! $acc& um_tau, vm_tau, residual_u, RU, res_u, ui, epsilon1, b_art, Re, artDivU, v_i, var,&
! $acc& vi, RV, res_v,t, v_i, Re, artDivV, var, residual_v, varM)


     CALL lacointerno(um_tau, um_n_tau, vm_tau, vm_n_tau, p, pn, T_tau, T_n_Tau)
       ! $acc data present(um_tau, um_n_tau, vm_tau, vm_n_tau, p, pn, T_tau, T_n_Tau) 
! $acc parallel      
 !       um_tau = um_n_tau
 !              vm_tau = vm_n_tau
 !              p      = pn
 !              T_tau = T_n_tau
! $acc end parallel          
  ! $acc end data

    !!--------- convergence criteria

!            IF (itc .NE. 1 .AND. error.LT.eps) then
!                EXIT
!            ENDIF
    
        ENDDO 
  
  !!------ END OF pseudo-time CALCULATION --------------------------------------------------- 
    CALL lacoexterno(um, vm, um_n_tau, vm_n_tau, T, T_n_tau)
 
!    IF (MOD(tr,n_tr).EQ.0) THEN
!        WRITE(*,*) '-----------------------------------------------------------'
!        WRITE(*,*) itc,'--Max Residual:',error
!        write(*,*) tr,'--Physical time:', time
!        WRITE(*,*) '-----------------------------------------------------------'

!        WRITE(*,*) '         dtau:',dtau   
!        WRITE(*,*) '         dt:',dt   
!        WRITE(*,*) '    Residual U:',residual_u
!        WRITE(*,*) '    Residual V:',residual_v
!        WRITE(*,*) '    Residual P:',residual_p    
!        WRITE(*,*) ' Artificial viscosity:',artMAX    
!        WRITE(*,*) ' Art Compressibility Par:',c2    
!        WRITE(*,*) '-----------------------------------------------------------'
!    ENDIF

        !CALL SYSTEM('rm data/error.dat')
    CALL iterat_convergence(itc, error)

       ! itc = 0
        !varM(1) = 0
       ! error = 100.d0
        !varM(2) = 100.d0
        !tr = tr + 1
        !varM(3) = varM(3) + 1

! $acc exit data copyout(um, um_n, vm, vm_n, RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p,&
! $acc& Z_tau, varZ, ym, xm, Pe, x, y, RZ, zi, T_tau, liga_poros, res_z, T, dtau, um_n_tau, vm_n_tau, T_n_tau, dx, dy,&
! $acc& um_tau, vm_tau, residual_u, RU, res_u, ui, epsilon1, b_art, Re, artDivU, v_i, var,&
! $acc& vi, RV, res_v,t, v_i, Re, artDivV, var, residual_v, varM)


!        open (550,file='data/time.dat')
!            
!             write (550,*) time
!
!        close(550)


    !-------- output preliminary results

!             IF (MOD(tr,n_tr).EQ.0) THEN
    
!                CALL comp_mean(u,v,um,vm)
                
!                CALL transient(u,v,p,T,Z,tr)

    !------ output data file ------------------------
!                CALL output(um,vm,u,v,p,Z,T,H,itc)
                
!            ENDIF

           ! tr = tr + 1
 
ENDDO


!$acc exit data copyout(um, um_n, vm, vm_n, RP, areau_e, areau_w, areav_n, areav_s, i, j, pi, p, dtau,c2, res_p, pn, residual_p,&
!$acc& Z_tau, varZ, ym, xm, Pe, x, y, RZ, zi, T_tau, liga_poros, res_z, T, dtau, um_n_tau, vm_n_tau, T_n_tau, dx, dy,&
!$acc& um_tau, vm_tau, residual_u, RU, res_u, ui, epsilon1, b_art, Re, artDivU, v_i, var,&
!$acc& vi, RV, res_v,t, v_i, Re, artDivV, var, residual_v, varM)



! $acc exit data copyout(um(1:imax+1,1:jmax),um_n(1:imax+1,1:jmax),res_u(1:imax+1,1:jmax),vm(1:imax,1:jmax+1),vm_n(1:imax,1:jmax+1),&
! $acc&res_v(1:imax,1:jmax+1),um_tau(1:imax+1,1:jmax),um_n_tau(1:imax+1,1:jmax),vm_tau(1:imax,1:jmax+1),vm_n_tau(1:imax,1:jmax+1),&
! $acc&u(1:imax,1:jmax),v(1:imax,1:jmax),P(1:imax,1:jmax),Pn(1:imax,1:jmax),H(1:imax,1:jmax),T(1:imax,1:jmax),Z(1:imax,1:jmax),&
! $acc&T_n_tau(1:imax,1:jmax),T_tau(1:imax,1:jmax),epsilon1(imax,jmax),areau_n(2:imax),areau_s(2:imax),areau_e(1:jmax),&
! $acc&areau_w(1:jmax),areav_n(1:imax),areav_s(1:imax),areav_e(2:jmax),areav_w(2:jmax),ui(1:imax+1,1:jmax),RU(1:imax+1,1:jmax),&
! $acc&RV(1:imax,1:jmax+1),vi(1:imax,1:jmax+1),&
! $acc&artDivU(imax,jmax),artDivV(imax,jmax),&
! $acc&x(1:imax),y(1:jmax),ym(1:jmax+1),xm(1:imax+1),dx, dy,&
! $acc&RP(1:imax,1:jmax),Pi(1:imax,1:jmax),res_p(1:imax,1:jmax),&
! $acc&RZ(1:imax,1:jmax),Zi(1:imax,1:jmax),res_z(2:imax-1,2:jmax-1),&
! $acc&Z_n_tau(1:imax,1:jmax),Z_tau(1:imax,1:jmax),liga_poros(imax,jmax), flag(imax,jmax), v_i,S, tsup, lf, q)

! $delete(um(1:imax+1,1:jmax),um_n(1:imax+1,1:jmax),res_u(1:imax+1,1:jmax),vm(1:imax,1:jmax+1),vm_n(1:imax,1:jmax+1),&
! $delete&res_v(1:imax,1:jmax+1),um_tau(1:imax+1,1:jmax),um_n_tau(1:imax+1,1:jmax),vm_tau(1:imax,1:jmax+1),vm_n_tau(1:imax,1:jmax+1),&
! $delete&u(1:imax,1:jmax),v(1:imax,1:jmax),P(1:imax,1:jmax),Pn(1:imax,1:jmax),H(1:imax,1:jmax),T(1:imax,1:jmax),Z(1:imax,1:jmax),&
! $delete&T_n_tau(1:imax,1:jmax),T_tau(1:imax,1:jmax),epsilon1(imax,jmax),areau_n(2:imax),areau_s(2:imax),areau_e(1:jmax),areau_w(1:jmax),&
! $delete&areav_n(1:imax),areav_s(1:imax),areav_e(2:jmax),areav_w(2:jmax),ui(1:imax+1,1:jmax),RU(1:imax+1,1:jmax),RV(1:imax,1:jmax+1),&
! $delete&vi(1:imax,1:jmax+1),&
! $delete&artDivU(imax,jmax),artDivV(imax,jmax),&
! $delete&x(1:imax),y(1:jmax),ym(1:jmax+1),xm(1:imax+1),dx,dy,&
! $delete&RP(1:imax,1:jmax),&Pi(1:imax,1:jmax),res_p(1:imax,1:jmax),&
! $delete&RZ(1:imax,1:jmax),Zi(1:imax,1:jmax),res_z(2:imax-1,2:jmax-1),&
! $delete&Z_n_tau(1:imax,1:jmax),&
! $delete&Z_tau(1:imax,1:jmax),liga_poros(imax,jmax), flag(imax,jmax), v_i, S, tsup, lf, q, Re, Pr, Pe)

!!------ END of physical CALCULATION -----------------------------------------------------

    open (550,file='data/time.dat')
            
        write (550,*) time

    close(550)

!!--------- calcular velocidades nos pontos medios
        CALL comp_mean(u,v,um,vm)
!!------------------------------------------------
            
        CALL transient(u,v,p,T,Z,itc)

!!------ output data file ------------------------
        CALL output(um,vm,u,v,p,Z,T,H,itc)

        clockTIME=ETIME(TARRAY)
		PRINT *, clockTime
		PRINT *, TARRAY
        WRITE(*,1009) clockTIME/60.d0!, 'horas' 
        1009  FORMAT(1X,'EXECUTION TIME = ',F15.4,' MINUTOS')

stop
END PROGRAM main

SUBROUTINE convergencia(varM, residual_u, residual_v, residual_p, itc)
    use comum
    use openacc
    !acc routine
    implicit none
    REAL(8), DIMENSION(10) :: varM
    REAL(8) :: residual_p, residual_u,residual_v, itc
    !$acc data present(itc, varM, residual_u, residual_v, residual_p)
        itc = itc+1
        varM(1) = varM(1) + 1
        varM(2) = MAX(residual_u,residual_v,residual_p)
    !$acc end data
    RETURN
END SUBROUTINE convergencia


SUBROUTINE lacoexterno(um, vm, um_n_tau, vm_n_tau, T, T_n_tau)
    use comum
    use openacc
    !acc routine
    implicit none
    REAL(8), DIMENSION(1:imax,1:jmax)     :: T_n_tau
    REAL(8), DIMENSION(1:imax+1,1:jmax) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) ::  um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) ::  vm_n_tau
    REAL(8), DIMENSION(1:imax,1:jmax)     :: T
    !REAL(8) :: error
    !REAL(8), DIMENSION(10) :: varM
    !INTEGER :: itc !, tr
!$acc data present(um, um_n_tau, vm, vm_n_tau, T, T_n_tau)
!$acc parallel
        um = um_n_tau 
        vm = vm_n_tau
        T = T_n_tau
!$acc end parallel
!$acc end data
        RETURN
END SUBROUTINE lacoexterno


SUBROUTINE lacointerno(um_tau, um_n_tau, vm_tau, vm_n_tau, p, pn, T_tau, T_n_Tau)
use comum
use openacc
!acc routine
    implicit none
REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, um_n_tau
REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, vm_n_tau
REAL(8), DIMENSION(1:imax,1:jmax)     :: P, Pn
REAL(8), DIMENSION(1:imax,1:jmax)     :: T_n_tau, T_tau

!$acc data present (um_tau, um_n_tau, vm_tau, vm_n_tau, p, pn, T_tau, T_n_tau)
!$acc parallel 
           um_tau = um_n_tau
            vm_tau = vm_n_tau
            p      = pn
            T_tau = T_n_tau
!$acc end parallel
!$acc end data
       RETURN
END SUBROUTINE lacointerno
