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
USE comum
IMPLICIT NONE

    INTEGER :: itc, tr, i, j 
    INTEGER*4 today(3), now(3)
    
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um, um_n, res_u
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm, vm_n, res_v
    
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, vm_n_tau

    REAL(8), DIMENSION(1:imax,1:jmax)     :: u, v, P, Pn, H, T,Z
    REAL(8), DIMENSION(1:imax,1:jmax)     :: T_n_tau, T_tau 
    REAL(8)                               :: InvFr2

    REAL(8) :: residual_p, residual_u,residual_v
    REAL(8) :: c2, error
    REAL(8) :: artMAX, yf, xf

    character(len=128) :: pwd
        REAL ETIME,clockTIME,TARRAY(2)
        clockTIME=ETIME(TARRAY)
        CALL idate(today)   ! today(1)=day, (2)=month, (3)=year
        CALL itime(now)     ! now(1)=hour, (2)=minute, (3)=second
        CALL get_environment_variable('PWD',pwd)
   
   
!------ input data -----------       
 
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

DO WHILE (time .LT. final_time)

     time = time + dt

    !!-------- pseudo-time calculation starts ---------------------------------------------------
        DO WHILE(itc.LT.itc_max)
          
    
    !!-------- Solve Momentum Equation with QUICK Scheme

            CALL solve_U(um,vm,um_n,um_tau,vm_tau,um_n_tau,pn,residual_u)
            
            CALL solve_V(um,vm,vm_n,um_tau,vm_tau,vm_n_tau,pn,T,InvFr2,residual_v)
    
    !!-------- Solve Continuity Equation 

            CALL solve_P(c2,p,um_n_tau,vm_n_tau,pn,residual_p)

    !-------- Solve Energy Equation ------------------------------------------------

            CALL solve_Z(um_n_tau,vm_n_tau,T,T_n_tau,T_tau)
            
    
       !     CALL solve_H(H,T_n_tau)
    
    
        !    CALL comp_T(Z_n_tau,H,T)
           
    !!-------- check convergence

            CALL convergence(itc,c2,error,residual_p,residual_u,residual_v)

            um_tau = um_n_tau
            vm_tau = vm_n_tau
            p      = pn
            T_tau = T_n_tau

    !!--------- convergence criteria
    
            IF (itc .NE. 1 .AND. error.LT.eps) then
            
                EXIT
            
            ENDIF
    
        ENDDO 
    !!------ END OF pseudo-time CALCULATION ---------------------------------------------------

        um = um_n_tau 
        vm = vm_n_tau
        T = T_n_tau

    IF (MOD(tr,n_tr).EQ.0) THEN
        WRITE(*,*) '-----------------------------------------------------------'
        WRITE(*,*) itc,'--Max Residual:',error
        write(*,*) tr,'--Physical time:', time
        WRITE(*,*) '-----------------------------------------------------------'

        WRITE(*,*) '         dtau:',dtau   
        WRITE(*,*) '         dt:',dt   
        WRITE(*,*) '    Residual U:',residual_u
        WRITE(*,*) '    Residual V:',residual_v
        WRITE(*,*) '    Residual P:',residual_p    
        WRITE(*,*) ' Artificial viscosity:',artMAX    
        WRITE(*,*) ' Art Compressibility Par:',c2    
        WRITE(*,*) '-----------------------------------------------------------'
    ENDIF

        !CALL SYSTEM('rm data/error.dat')
        itc = 0
        error = 100.d0 

!        open (550,file='data/time.dat')
!            
!             write (550,*) time
!
!        close(550)


    !-------- output preliminary results

             IF (MOD(tr,n_tr).EQ.0) THEN
    
                CALL comp_mean(u,v,um,vm)
                
                CALL transient(u,v,p,T,Z,tr)

    !------ output data file ------------------------
                CALL output(um,vm,u,v,p,Z,T,H,itc)
                
            ENDIF

            tr = tr + 1
 
ENDDO
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
        WRITE(*,1009) clockTIME/60.d0!, 'horas' 
        1009  FORMAT(1X,'EXECUTION TIME = ',F15.4,' MINUTOS')

stop
END PROGRAM main


