MODULE comum

        NAMELIST /iterations/ itc_max, nc, n_tr, n_out, n_vort,&
                              beta, b_art, dtau_f,final_time, eps, eps_mass, start_mode
        NAMELIST /ref/ nu, YF_b, YO_oo, Ts, nToo

        INTEGER  :: itc_max        !numero de iterações
        
        !frequencia dos outputs:
        INTEGER  :: nc             !erros
        INTEGER  :: n_tr           !plota parte transiente
        INTEGER  :: n_out          !salva resultados preliminares
        INTEGER  :: n_vort         !salva dados do vortice
        INTEGER  :: restart_mode   !tipo de start, se eh CI ou solucao anterior
                
        REAL(8)  :: dtau, dt, dt_p, time, pseudo_time, final_time        !passo de tempo
        REAL(8)  :: eps , eps_mass  !criterio de convergencia  
        REAL(8)  :: beta          !parâmetro de compressi        
        REAL(8)  :: b_art         !coeficiente da dissipação artificial
        REAL(8)  :: dtau_f          !fator de correcao para calc de dt        
        REAL(8), PARAMETER  :: porosidade=0.2d0
        REAL(8), PARAMETER  :: Darcy_number = 1.0d-2
        REAL(8), PARAMETER  :: Cf = 1.d0
    !    REAL(8)  :: Liga_poros
        REAL(8), PARAMETER  :: Temp_cylinder = 5.d0
   !     REAL(8)  :: Darcy_term
        !parametros de refinamento P e Q
        !colocar P=1 para a malha uniforme
        REAL(8), PARAMETER :: Px_grid = 1.6d0, Py_grid = 1.6d0, Q_grid = 1.8d0 
       
        !tamanho do dominio 
        REAL(8), PARAMETER :: Lhori = 12.    !largura total
        
        REAL(8), PARAMETER :: y_up = 10      !altura em y+
       
        REAL(8), PARAMETER :: y_down = 5.   !altura em y-

        REAL(8), PARAMETER :: Hvert = y_up + y_down !altura toral

        INTEGER, PARAMETER :: imax  = 51 !numero de pontos da malha em x

        REAL(8), PARAMETER :: dx_c = Lhori / (imax-1) !dita o tamanho de dy e dx

!        INTEGER, PARAMETER :: imax  = int_points + 1! (Lhori / dx_c) + 1!numero de pontos da malha em x

        INTEGER, PARAMETER :: jmax  = (Hvert / dx_c) + 1!numero de pontos da malha em y

        REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
        REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada

        REAL(8)            :: vol_u(2:imax,1:jmax)       !volume de controle de u
        REAL(8)            :: vol_v(1:imax,2:jmax)       !volume de controle de v
        REAL(8)            :: vol_p(1:imax,1:jmax)       !volume de controle de p

        REAL(8)            :: areau_n(2:imax)            !area n de u
        REAL(8)            :: areau_s(2:imax)            !area s de u
        REAL(8)            :: areau_e(1:jmax)            !area e de u
        REAL(8)            :: areau_w(1:jmax)            !area w de u

        REAL(8)            :: areav_n(1:imax)            !area n de v
        REAL(8)            :: areav_s(1:imax)            !area s de v
        REAL(8)            :: areav_e(2:jmax)            !area e de v
        REAL(8)            :: areav_w(2:jmax)            !area w de v
        
        REAL(8) 		:: dx(2:imax), dy(2:jmax)
 	     REAL(8) 		:: epsilon1(imax,jmax)
        REAL(8)            	:: liga_poros(imax,jmax)   
        REAL(8), PARAMETER :: rad1 = 1.d0 !raio do cilindro        

    
        ! flags for obstacle interior, boundary, fluid cells and close to boundary
        INTEGER, PARAMETER :: C_I = 2, C_B = 1, C_F = 0, C_BS = 3 
          
        INTEGER                :: flag(imax,jmax)

        REAL(8)                :: g = 9.80665d0         !gravitational constant [m/s^2]
        REAL(8)                :: ao  = 1.d-3           !initial radius [m]
        REAL(8)                :: v_idim = 1.55d-1       !velocidade de injecao [m/s]
        REAL(8)                :: q_f = 0.413995523d0
        REAL(8)                :: L_c
        REAL(8)                :: v_c
        REAL(8)                :: v_i
        REAL(8)                :: Fr
        REAL(8)                :: S
        REAL(8)                :: Lf = 1.d0
        REAL(8)                :: Lo = 1.d0
        REAL(8)                :: nu    ! for Methane , 3.51d0 for n-Heptane
        REAL(8)                :: YF_b 
        REAL(8)                :: YO_oo 
        REAL(8)                :: Ts    ! Tb  = boiling temperature [k]
        REAL(8)                :: Too   ! Too = dimen ambient temp [k]
        REAL(8)                :: nToo
        REAL(8)                :: Tsup
        REAL(8)                :: Tinf

!        REAL(8)                :: q_dim = 47904.1916d3        ! q    = combustion heat [J/kg] for n-Heptane C7H16
        REAL(8)                :: q_dim = 50.15d6        ! q    = combustion heat [J/kg] for Methane CH4
        REAL(8)                :: q  
        REAL(8)                :: cp_tot    != 1937.3540735d0! [J/kgK]
        REAL(8)                :: rho_tot   != 1.1950251341d0  ! [kg/m^3]
        REAL(8)                :: k_tot     != 7.3404587d-2    ! [W/mK]
        REAL(8)                :: nu_tot    != 9.30492d-5      ! [m^2/s]
        REAL(8)                :: alpha_tot                   ! [m^2/s]
        REAL(8)                :: Re ,Pr, Pe                 

end MODULE comum
