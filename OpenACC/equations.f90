subroutine RESU(um_tau,vm_tau,pn,RU)
    use comum
    ! $acc routine
    implicit none
    integer :: i, j   
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: Pn
    !REAL(8), DIMENSION(1:imax+1,1:jmax+1) ::epsilon1
    real(8) :: fw,fe,fs,fn, df, aw,aww,ae,aee,as,ass,an,ann, ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: afw, afe, afn, afs
    real(8) :: q_art !artDivU(imax,jmax)
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
    !$acc data present(um_tau,vm_tau,Pn,RU,epsilon1,x,y,xm,ym,artDivU,liga_poros,areau_s,areau_w,areau_n,areau_e)

    ! $omp parallel do private(q_art, i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dudxdx,dxdvdy,q_art,u_W,u_WW,u_E,u_EE,u_N,u_NN,u_S,u_SS,v_P,u_P,ap,aw,aww,ae,aee,an,ann,as,ass)

    !$acc parallel loop collapse(2) &
    !$acc&private(fw,dw,ds,v_p,u_ww,u_nn,u_n,u_ee,u_s,u_p,u_e,q_art,u_w,u_ss,ann,an,as,ae,aw,ap,ass,afs,afn,afe,aee,afw,dudxdx,aww,dn,df,de,fs,dxdvdy,fn,fe,i,j)
    DO j=3,jmax-2
      DO i=3,imax-1
        fn = 0.5d0 * ( vm_tau(i  ,j+1) + vm_tau(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        fs = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        fe = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        fw = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )

        df = fe - fw + fn - fs

        Dn = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        Ds = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        De = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        Dw = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))

        !CALL quick(ap,aw,aww,ae,aee,an,ann,as,ass,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
        if(fw.GT.0.0d0) then
                afw = 1.0d0
        elseif(fw.LT.0.0d0) then
                afw = 0.0d0
        endif

        if(fe.GT.0.0d0) then
                afe = 1.0d0
        elseif(fe.LT.0.0d0) then
                afe = 0.0d0
        endif

        if(fn.GT.0.0d0) then
                afn = 1.0d0
        elseif(fn.LT.0.0d0) then
                afn = 0.0d0
        endif

        if(fs.GT.0.0d0) then
                afs = 1.0d0
        elseif(fs.LT.0.0d0) then
                afs = 0.0d0
        endif

        aw = Dw + 0.75d0  * afw * fw &
                + 0.125d0 * afe * fe &
                + 0.375d0 * ( 1.0d0 - afw ) * fw

        ae = De - 0.375d0 * afe * fe &
                - 0.75d0  * ( 1.0d0 - afe ) * fe &
                - 0.125d0 * ( 1.0d0 - afw ) * fw

        as = Ds + 0.75d0  * afs * fs &
                + 0.125d0 * afn * fn &
                + 0.375d0 * ( 1.0d0 - afs ) * fs

        an = Dn - 0.375d0 * afn * fn &
                - 0.75d0  * ( 1.0d0 - afn ) * fn &
                - 0.125d0 * ( 1.0d0 - afs ) * fs

        aww = -0.125d0 *           afw   * fw
        aee =  0.125d0 * ( 1.0d0 - afe ) * fe
        ass = -0.125d0 *           afs   * fs
        ann =  0.125d0 * ( 1.0d0 - afn ) * fn

        ap = aw + ae + as + an + aww + aee + ass + ann + df
        !end Quick

        u_W  = um_tau(i-1,j  )
        u_WW = um_tau(i-2,j  )
        u_E  = um_tau(i+1,j  )
        u_EE = um_tau(i+2,j  )
        u_S  = um_tau(i  ,j-1)
        u_SS = um_tau(i  ,j-2)
        u_N  = um_tau(i  ,j+1)
        u_NN = um_tau(i  ,j+2)
        u_P = um_tau(i  ,j  )
        v_P = vm_tau(i  ,j  )

        dudxdx =  areau_e(j) * ( u_E - u_P ) / (xm(i+1)-xm(i  )) &
                 -areau_w(j) * ( u_P - u_W ) / (xm(i  )-xm(i-1)) 

        dxdvdy =  areau_e(j) * (vm_tau(i  ,j+1) - vm_tau(i  ,j  )) / (ym(j+1)-ym(j)) &
                 -areau_w(j) * (vm_tau(i-1,j+1) - vm_tau(i-1,j  )) / (ym(j+1)-ym(j))

        artDivU(i,j) = - b_art *(dudxdx + dxdvdy)

        !bulk artificial viscosity term from Ramshaw(1990)
        q_art =epsilon1(i,j)*  ( pn(i,j)-pn(i-1,j) )/ (x(i)-x(i-1)) + artDivU(i,j)

        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * u_P &
                          + aww * u_WW + aw * u_W &
                          + aee * u_EE + ae * u_E &
                          + ass * u_SS + as * u_S &
                          + ann * u_NN + an * u_N)&
                          - q_art &
            -epsilon1(i,j)*(u_P/(Re*Darcy_number) + &
            CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* u_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
      ENDDO
    ENDDO
!$acc end data
RETURN
END SUBROUTINE RESU

!#######################################################################################


!#######################################################################################
subroutine RESU1(ui,vm_tau,pn,RU)
    use comum
    ! $acc routine
    implicit none
    integer :: i, j   
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: ui, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: Pn
    !REAL(8), DIMENSION(1:imax+1,1:jmax+1) ::epsilon1
    real(8) :: fw,fe,fs,fn, df, aw,aww,ae,aee,as,ass,an,ann, ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: afw, afe, afn, afs
    real(8) :: q_art !artDivU(imax,jmax)
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
    !$acc data present(ui,vm_tau,Pn,RU,epsilon1,x,y,xm,ym,artDivU,liga_poros,areau_s,areau_w,areau_n,areau_e)

    ! $omp parallel do private(q_art, i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dudxdx,dxdvdy,q_art,u_W,u_WW,u_E,u_EE,u_N,u_NN,u_S,u_SS,v_P,u_P,ap,aw,aww,ae,aee,an,ann,as,ass)

    !$acc parallel loop collapse(2) &
    !$acc&private(fw,dw,ds,v_p,u_ww,u_nn,u_n,u_ee,u_s,u_p,u_e,q_art,u_w,u_ss,ann,an,as,ae,aw,ap,ass,afs,afn,afe,aee,afw,dudxdx,aww,dn,df,de,fs,dxdvdy,fn,fe,i,j)
    DO j=3,jmax-2
      DO i=3,imax-1
        fn = 0.5d0 * ( vm_tau(i  ,j+1) + vm_tau(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        fs = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        fe = 0.5d0 * ( ui(i+1,j  ) + ui(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        fw = 0.5d0 * ( ui(i  ,j  ) + ui(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )

        df = fe - fw + fn - fs

        Dn = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        Ds = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        De = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        Dw = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))

        !CALL quick(ap,aw,aww,ae,aee,an,ann,as,ass,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
        if(fw.GT.0.0d0) then
                afw = 1.0d0
        elseif(fw.LT.0.0d0) then
                afw = 0.0d0
        endif

        if(fe.GT.0.0d0) then
                afe = 1.0d0
        elseif(fe.LT.0.0d0) then
                afe = 0.0d0
        endif

        if(fn.GT.0.0d0) then
                afn = 1.0d0
        elseif(fn.LT.0.0d0) then
                afn = 0.0d0
        endif

        if(fs.GT.0.0d0) then
                afs = 1.0d0
        elseif(fs.LT.0.0d0) then
                afs = 0.0d0
        endif

        aw = Dw + 0.75d0  * afw * fw &
                + 0.125d0 * afe * fe &
                + 0.375d0 * ( 1.0d0 - afw ) * fw

        ae = De - 0.375d0 * afe * fe &
                - 0.75d0  * ( 1.0d0 - afe ) * fe &
                - 0.125d0 * ( 1.0d0 - afw ) * fw

        as = Ds + 0.75d0  * afs * fs &
                + 0.125d0 * afn * fn &
                + 0.375d0 * ( 1.0d0 - afs ) * fs

        an = Dn - 0.375d0 * afn * fn &
                - 0.75d0  * ( 1.0d0 - afn ) * fn &
                - 0.125d0 * ( 1.0d0 - afs ) * fs

        aww = -0.125d0 *           afw   * fw
        aee =  0.125d0 * ( 1.0d0 - afe ) * fe
        ass = -0.125d0 *           afs   * fs
        ann =  0.125d0 * ( 1.0d0 - afn ) * fn

        ap = aw + ae + as + an + aww + aee + ass + ann + df
        !end Quick

        u_W  = ui(i-1,j  )
        u_WW = ui(i-2,j  )
        u_E  = ui(i+1,j  )
        u_EE = ui(i+2,j  )
        u_S  = ui(i  ,j-1)
        u_SS = ui(i  ,j-2)
        u_N  = ui(i  ,j+1)
        u_NN = ui(i  ,j+2)
        u_P = ui(i  ,j  )
        v_P = vm_tau(i  ,j  )

        dudxdx =  areau_e(j) * ( u_E - u_P ) / (xm(i+1)-xm(i  )) &
                 -areau_w(j) * ( u_P - u_W ) / (xm(i  )-xm(i-1)) 

        dxdvdy =  areau_e(j) * (vm_tau(i  ,j+1) - vm_tau(i  ,j  )) / (ym(j+1)-ym(j)) &
                 -areau_w(j) * (vm_tau(i-1,j+1) - vm_tau(i-1,j  )) / (ym(j+1)-ym(j))

        artDivU(i,j) = - b_art *(dudxdx + dxdvdy)   

        !bulk artificial viscosity term from Ramshaw(1990)
        q_art =epsilon1(i,j)*  ( pn(i,j)-pn(i-1,j) )/ (x(i)-x(i-1)) + artDivU(i,j)

        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * u_P &
                          + aww * u_WW + aw * u_W &
                          + aee * u_EE + ae * u_E &
                          + ass * u_SS + as * u_S &
                          + ann * u_NN + an * u_N)&
                          - q_art &
            -epsilon1(i,j)*(u_P/(Re*Darcy_number) + &
            CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* u_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
      ENDDO
    ENDDO
!$acc end data
RETURN
END SUBROUTINE RESU1

!#######################################################################################
subroutine upwind_U(um_tau,vm_tau,pn,RU,ii,iii,jj,jjj)
    use comum
    ! $acc routine
    implicit none
    integer :: i, j
    integer :: ii, iii,jj, jjj    
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: Pn
    !REAL(8), DIMENSION(1:imax+1,1:jmax+1) ::epsilon1
    real(8) :: fw,fe,fs,fn, df, aw,aww,ae,aee,as,ass,an,ann, ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: q_art
    real(8) :: dudxdx, dxdvdy
    !$acc data present(ru(:,:),liga_poros(:,:),pn(:,:),areau_e(:),areau_n(:),areau_s(:),areau_w(:),vm_tau(:,:),xm(:),x(:),ym(:),y(:),epsilon1(:,:),ui(:,:),um_tau)

    !!!!!!!!!!!!!!!!!!!!!compute x-direction velocity component un!!!!!!!!!!!!!!!!!!!!!
    ! $omp parallel do private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dudxdx,dxdvdy,q_art,u_W,u_E,u_N,u_S,u_P,V_P,ap,aw,ae,an,as)

    !$acc parallel loop collapse(2) private(fw,dw,ds,v_p,u_w,u_n,u_e,u_p,q_art,u_s,an,ae,ap,as,dudxdx,aw,dn,df,de,fs,dxdvdy,fn,fe,i,j)
    do j=jj,jjj
      do i=ii,iii
        fn = 0.5d0 * ( vm_tau(i  ,j+1) + vm_tau(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        fs = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        fe = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        fw = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )

        df = fe - fw + fn - fs

        Dn = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        Ds = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        De = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        Dw = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))

        !CALL upwind(ap,aw,ae,an,as,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
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
            
        dxdvdy =  areau_e(j) * (vm_tau(i  ,j+1) - vm_tau(i  ,j  )) / (ym(j+1)-ym(j)) &
                 -areau_w(j) * (vm_tau(i-1,j+1) - vm_tau(i-1,j  )) / (ym(j+1)-ym(j))

        !!bulk artificial viscosity term from Ramshaw(1990)
        q_art = epsilon1(i,j)* ( pn(i,j)-pn(i-1,j) )/ (x(i)-x(i-1)) - b_art * (dudxdx + dxdvdy)

        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * u_P &
                          + aw * u_W + ae * u_E  &
                          + as * u_S + an * u_N) &
                          - q_art &
            - epsilon1(i,j)* (  u_P/(Re*Darcy_number) +&
            CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* u_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
    enddo
  enddo
!$acc end data
RETURN
END SUBROUTINE upwind_U


!#######################################################################################
subroutine upwind_U1(ui,vm_tau,pn,RU,ii,iii,jj,jjj)
    use comum
    ! $acc routine
    implicit none
    integer :: i, j
    integer :: ii, iii,jj, jjj    
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: ui, RU
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: Pn
    !REAL(8), DIMENSION(1:imax+1,1:jmax+1) ::epsilon1
    real(8) :: fw,fe,fs,fn, df, aw,aww,ae,aee,as,ass,an,ann, ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: q_art
    real(8) :: dudxdx, dxdvdy
    !$acc data present(ru(:,:),liga_poros(:,:),pn(:,:),areau_e(:),areau_n(:),areau_s(:),areau_w(:),vm_tau(:,:),xm(:),x(:),ym(:),y(:),epsilon1(:,:),ui(:,:))

    !!!!!!!!!!!!!!!!!!!!!compute x-direction velocity component un!!!!!!!!!!!!!!!!!!!!!
    ! $omp parallel do private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dudxdx,dxdvdy,q_art,u_W,u_E,u_N,u_S,u_P,V_P,ap,aw,ae,an,as)

    !$acc parallel loop collapse(2) private(fw,dw,ds,v_p,u_w,u_n,u_e,u_p,q_art,u_s,an,ae,ap,as,dudxdx,aw,dn,df,de,fs,dxdvdy,fn,fe,i,j)
    do j=jj,jjj
      do i=ii,iii
        fn = 0.5d0 * ( vm_tau(i  ,j+1) + vm_tau(i-1,j+1) ) * areau_n(i) / (epsilon1(i,j) )
        fs = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i-1,j  ) ) * areau_s(i) / (epsilon1(i,j) )
        fe = 0.5d0 * ( ui(i+1,j  ) + ui(i  ,j  ) ) * areau_e(j) / (epsilon1(i,j) )
        fw = 0.5d0 * ( ui(i  ,j  ) + ui(i-1,j  ) ) * areau_w(j) / (epsilon1(i,j) )

        df = fe - fw + fn - fs

        Dn = (epsilon1(i,j)/Re) * areau_n(i) / (y(j+1)-y(j  ))
        Ds = (epsilon1(i,j)/Re) * areau_s(i) / (y(j  )-y(j-1))
        De = (epsilon1(i,j)/Re) * areau_e(j) / (xm(i+1)-xm(i  ))
        Dw = (epsilon1(i,j)/Re) * areau_w(j) / (xm(i  )-xm(i-1))

        !CALL upwind(ap,aw,ae,an,as,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
        aw = Dw + MAX(fw , 0.0d0)
        as = Ds + MAX(fs , 0.0d0)
    
        ae = De + MAX(0.0d0 , -fe)
        an = Dn + MAX(0.0d0 , -fn) 

        ap = aw + ae + as + an + df
    
        u_W = ui(i-1,j  )
        u_E = ui(i+1,j  )
        u_S = ui(i  ,j-1)
        u_N = ui(i  ,j+1)
        u_P = ui(i  ,j  )
        v_P = vm_tau(i  ,j  )
        dudxdx =  areau_e(j) * ( u_E - u_P ) / (xm(i+1)-xm(i  )) &
                 -areau_w(j) * ( u_P - u_W ) / (xm(i  )-xm(i-1)) 
            
        dxdvdy =  areau_e(j) * (vm_tau(i  ,j+1) - vm_tau(i  ,j  )) / (ym(j+1)-ym(j)) &
                 -areau_w(j) * (vm_tau(i-1,j+1) - vm_tau(i-1,j  )) / (ym(j+1)-ym(j))

        !!bulk artificial viscosity term from Ramshaw(1990)
        q_art = epsilon1(i,j)* ( pn(i,j)-pn(i-1,j) )/ (x(i)-x(i-1)) - b_art * (dudxdx + dxdvdy)

        RU(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * u_P &
                          + aw * u_W + ae * u_E  &
                          + as * u_S + an * u_N) &
                          - q_art &
            - epsilon1(i,j)* (  u_P/(Re*Darcy_number) +&
            CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* u_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
    enddo
  enddo
!$acc end data
RETURN
END SUBROUTINE upwind_U1

!#######################################################################################
subroutine solve_U(um,vm,um_n,um_tau,vm_tau,um_n_tau,pn,residual_u)
!subroutine solve_U(um,vm,um_n,um_tau,vm_tau,um_n_tau,p,residual_u,res_u,RU,UI,areau_e,areau_w,areau_n,areau_s,ym,xm,Re,x,y,liga_poros,dtau,epsilon1,artDivU)
    use comum
    !use omp_lib
    ! $acc routine
    ! $acc routine(RESU)
    ! $acc routine(RESU1)
    ! $acc routine(upwind_U)
    ! $acc routine(upwind_U1)
    !$acc routine(bcUV)
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um,um_n!,ui,RU
    !REAL(8), DIMENSION(3:imax-1,2:jmax-1) :: res_u
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau, um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax  ) :: pn
    REAL(8) :: residual_u
    !REAL(8)            :: artDivU(imax,jmax)
    !REAL(8)            :: areau_n(2:imax)            !area n de u
    !REAL(8)            :: areau_s(2:imax)            !area s de u
    !REAL(8)            :: areau_e(1:jmax)            !area e de u
    !REAL(8)            :: areau_w(1:jmax)            !area w de u
    !REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    !REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    !REAL(8)            :: Re,dtau
    !REAL(8)            :: liga_poros(imax,jmax)
    !REAL(8), DIMENSION (1:imax,1:jmax) :: epsilon1

    !$acc data present(res_u,um_n_tau,epsilon1,liga_poros,areau_e,areau_n,areau_w,um_tau,x,vm_tau,xm,y,ym,areau_s,dtau,pn,ui,ru,re,artdivu,um,vm)
    !RU = 0.d0
    !res_u = 0.d0

    write (*,*) vm

    CALL RESU(um_tau,vm_tau,pn,RU)

    ! $omp parallel do private(i) 
    !DO i=2,imax
      CALL upwind_U(um_tau,vm_tau,pn,RU,2,imax,2,2)
    !ENDDO

    ! $omp parallel do private(i) 
    !DO i=2,imax
      CALL upwind_U(um_tau,vm_tau,pn,RU,2,imax,jmax-1,jmax-1)
    !ENDDO

    ! $omp parallel do private(j) 
    !DO j=2,jmax-1
      CALL upwind_U(um_tau,vm_tau,pn,RU,2,2,2,jmax-1)
    !ENDDO

    ! $omp parallel do private(j) 
    !DO j=2,jmax-1
      CALL upwind_U(um_tau,vm_tau,pn,RU,imax,imax,2,jmax-1)
    !ENDDO

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j)
    DO j=2,jmax-1
      DO i=3,imax-1
        res_u(i,j) =( (um(i,j)-um_tau(i,j)) +  RU(i,j)*dt) * dtau 
        ui(i,j) =  ( um_tau(i,j) + res_u(i,j))
      ENDDO
    ENDDO

    call bcUV1(ui,vm_tau)

    CALL RESU1(ui,vm_tau,pn,RU)

    ! $omp parallel do private(i) 
    !DO i=2,imax
      CALL upwind_U1(ui,vm_tau,pn,RU,2,imax,2,2)
    !ENDDO

    ! $omp parallel do private(i) 
    !DO i=2,imax
      CALL upwind_U1(ui,vm_tau,pn,RU,2,imax,jmax-1,jmax-1)
    !ENDDO

    ! $omp parallel do private(j) 
    !DO j=2,jmax-1
      CALL upwind_U1(ui,vm_tau,pn,RU,2,2,2,jmax-1)
    !ENDDO

    ! $omp parallel do private(j) 
    !DO j=2,jmax-1
      CALL upwind_U1(ui,vm_tau,pn,RU,imax,imax,2,jmax-1)
    !ENDDO

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j)
    DO j=2,jmax-1
      DO i=3,imax-1
        res_u(i,j) =( (um(i,j)-um_tau(i,j)) + RU(i,j)*dt) * dtau 
        ui(i,j) = (0.75d0 * um_tau(i,j) + 0.25d0 * &
                  ( ui(i,j) + res_u(i,j) ) )
      ENDDO
    ENDDO

    call bcUV1(ui,vm_tau)

    CALL RESU1(ui,vm_tau,pn,RU)

    ! $omp parallel do private(i) 
    !DO i=2,imax
      CALL upwind_U1(ui,vm_tau,pn,RU,2,imax,2,2)
    !ENDDO

    ! $omp parallel do private(i) 
    !DO i=2,imax
      CALL upwind_U1(ui,vm_tau,pn,RU,2,imax,jmax-1,jmax-1)
    !ENDDO

    ! $omp parallel do private(j) 
    !DO j=2,jmax-1
      CALL upwind_U1(ui,vm_tau,pn,RU,2,2,2,jmax-1)
    !ENDDO

    ! $omp parallel do private(j) 
    !DO j=2,jmax-1
      CALL upwind_U1(ui,vm_tau,pn,RU,imax,imax,2,jmax-1)
    !ENDDO

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j)
    DO j=2,jmax-1
      DO i=3,imax-1
        res_u(i,j) =( (um(i,j)-um_tau(i,j)) +  RU(i,j)*dt) * dtau 
        um_n_tau(i,j) = (1.0d0 / 3.0d0 * um_tau(i,j) + 2.0d0 / 3.0d0 * &
                        ( ui(i,j) + res_u(i,j) ) )
      ENDDO
    ENDDO

    call bcUV2(um_n_tau,vm_tau)
    
    !$acc parallel
    residual_u =  MAXVAL(ABS(res_u)) 
    !$acc end parallel

!    open (550,file='data/resu.dat')        
!         do i = 1,imax
!             do j = 1,jmax
!             write (550,*) X(i) , Y(j) , res_u(i,j)
!             enddo
!         enddo
!    close(550)
! $acc end kernels 
!$acc end data
RETURN
END SUBROUTINE solve_U

!#######################################################################################
subroutine RESV(um_tau,vm_tau,pn,T,InvFr2,RV)
    use comum
    ! $acc routine
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: Pn,T
    !REAL(8), DIMENSION(1:imax,1:jmax) :: epsilon1
    real(8) :: InvFr2
    real(8) :: fw,fe,fs,fn, df, aw,aww,ae,aee,as,ass,an,ann, ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: afw, afe, afn, afs
    real(8) :: q_art!, artDivV(imax,jmax)! artDivU(imax,jmax)
    real(8) :: dvdydy, dydudx
    !REAL(8)            :: areav_n(1:imax)            !area n de v
    !REAL(8)            :: areav_s(1:imax)            !area s de v
    !REAL(8)            :: areav_e(2:jmax)            !area e de v
    !REAL(8)            :: areav_w(2:jmax)            !area w de v
    !REAL(8)            :: Re
    !REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    !REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    !REAL(8)            :: b_art         !coeficiente da dissipação artificial
    !REAL(8)            :: liga_poros(imax,jmax)

    !$acc data present(um_tau,vm_tau,res_v,pn,t,liga_poros,rv,areav_s,areav_n,areav_e,areav_w,artdivv,ym,y,x,xm,epsilon1)
    ! $omp parallel do private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dvdydy,dydudx,q_art,v_W,v_WW,v_E,v_EE,v_N,v_NN,v_S,v_SS,v_P,u_P,ap,aw,aww,ae,aee,an,ann,as,ass)

    !$acc parallel loop collapse(2) private(dydudx,dw,ds,q_art,v_nn,v_n,v_ee,v_s,v_p,v_e,u_p,v_ww,v_w,v_ss,ann,an,as,ae,aw,ap,ass,afs,afn,afe,aee,afw,dvdydy,aww,dn,df,de,fw,fs,fn,fe)
    DO j=3,jmax-1
      DO i=3,imax-2
        fn = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        fs = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        fe = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        fw = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)

        df = fe - fw + fn - fs

        Dn = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))

        !CALL quick(ap,aw,aww,ae,aee,an,ann,as,ass,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
        if(fw.GT.0.0d0) then
                afw = 1.0d0
        elseif(fw.LT.0.0d0) then
                afw = 0.0d0
        endif

        if(fe.GT.0.0d0) then
                afe = 1.0d0
        elseif(fe.LT.0.0d0) then
                afe = 0.0d0
        endif

        if(fn.GT.0.0d0) then
                afn = 1.0d0
        elseif(fn.LT.0.0d0) then
                afn = 0.0d0
        endif

        if(fs.GT.0.0d0) then
                afs = 1.0d0
        elseif(fs.LT.0.0d0) then
                afs = 0.0d0
        endif

        aw = Dw + 0.75d0  * afw * fw &
                + 0.125d0 * afe * fe &
                + 0.375d0 * ( 1.0d0 - afw ) * fw

        ae = De - 0.375d0* afe * fe &
                - 0.75d0  * ( 1.0d0 - afe ) * fe &
                - 0.125d0 * ( 1.0d0 - afw ) * fw

        as = Ds + 0.75d0  * afs * fs &
                + 0.125d0 * afn * fn &
                + 0.375d0 * ( 1.0d0 - afs ) * fs

        an = Dn - 0.375d0* afn * fn &
                - 0.75d0  * ( 1.0d0 - afn ) * fn &
                - 0.125d0 * ( 1.0d0 - afs ) * fs

        aww = -0.125d0 *           afw   * fw
        aee =  0.125d0 * ( 1.0d0 - afe ) * fe
        ass = -0.125d0 *           afs   * fs
        ann =  0.125d0 * ( 1.0d0 - afn ) * fn

        ap = aw + ae + as + an + aww + aee + ass + ann + df
        !end Quick

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
            
        dydudx =  areav_n(i) * (um_tau(i+1,j  )-um_tau(i,j  )) / (xm(i+1)-xm(i)) &
                 -areav_s(i) * (um_tau(i+1,j-1)-um_tau(i,j-1)) / (xm(i+1)-xm(i))

        artDivV(i,j) = - b_art *(dydudx + dvdydy)            

        !bulk artificial viscosity term from Ramshaw(1990)
        q_art = epsilon1(i,j)* ( pn(i,j)-pn(i,j-1) )/ (y(j)-y(j-1)) + artDivV(i,j)
 
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
!$acc end data
RETURN
END SUBROUTINE RESV


!#######################################################################################
subroutine RESV1(um_tau,vi,pn,T,InvFr2,RV)
    use comum
    ! $acc routine
    implicit none
    integer :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vi, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: Pn,T
    !REAL(8), DIMENSION(1:imax,1:jmax) :: epsilon1
    real(8) :: InvFr2
    real(8) :: fw,fe,fs,fn, df, aw,aww,ae,aee,as,ass,an,ann, ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: afw, afe, afn, afs
    real(8) :: q_art!, artDivV(imax,jmax)! artDivU(imax,jmax)
    real(8) :: dvdydy, dydudx
    !REAL(8)            :: areav_n(1:imax)            !area n de v
    !REAL(8)            :: areav_s(1:imax)            !area s de v
    !REAL(8)            :: areav_e(2:jmax)            !area e de v
    !REAL(8)            :: areav_w(2:jmax)            !area w de v
    !REAL(8)            :: Re
    !REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    !REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    !REAL(8)            :: b_art         !coeficiente da dissipação artificial
    !REAL(8)            :: liga_poros(imax,jmax)

    !$acc data present(um_tau,res_v,pn,t,vi,liga_poros,rv,areav_s,areav_n,areav_e,areav_w,artdivv,ym,y,x,xm,epsilon1)
    ! $omp parallel do private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dvdydy,dydudx,q_art,v_W,v_WW,v_E,v_EE,v_N,v_NN,v_S,v_SS,v_P,u_P,ap,aw,aww,ae,aee,an,ann,as,ass)

    !$acc parallel loop collapse(2) private(dydudx,dw,ds,q_art,v_nn,v_n,v_ee,v_s,v_p,v_e,u_p,v_ww,v_w,v_ss,ann,an,as,ae,aw,ap,ass,afs,afn,afe,aee,afw,dvdydy,aww,dn,df,de,fw,fs,fn,fe)
    DO j=3,jmax-1
      DO i=3,imax-2
        fn = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j+1) ) * areav_n(i) / epsilon1(i,j)
        fs = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)
        fe = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)
        fw = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)

        df = fe - fw + fn - fs

        Dn = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))

        !CALL quick(ap,aw,aww,ae,aee,an,ann,as,ass,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
        if(fw.GT.0.0d0) then
                afw = 1.0d0
        elseif(fw.LT.0.0d0) then
                afw = 0.0d0
        endif

        if(fe.GT.0.0d0) then
                afe = 1.0d0
        elseif(fe.LT.0.0d0) then
                afe = 0.0d0
        endif

        if(fn.GT.0.0d0) then
                afn = 1.0d0
        elseif(fn.LT.0.0d0) then
                afn = 0.0d0
        endif

        if(fs.GT.0.0d0) then
                afs = 1.0d0
        elseif(fs.LT.0.0d0) then
                afs = 0.0d0
        endif

        aw = Dw + 0.75d0  * afw * fw &
                + 0.125d0 * afe * fe &
                + 0.375d0 * ( 1.0d0 - afw ) * fw

        ae = De - 0.375d0* afe * fe &
                - 0.75d0  * ( 1.0d0 - afe ) * fe &
                - 0.125d0 * ( 1.0d0 - afw ) * fw

        as = Ds + 0.75d0  * afs * fs &
                + 0.125d0 * afn * fn &
                + 0.375d0 * ( 1.0d0 - afs ) * fs

        an = Dn - 0.375d0* afn * fn &
                - 0.75d0  * ( 1.0d0 - afn ) * fn &
                - 0.125d0 * ( 1.0d0 - afs ) * fs

        aww = -0.125d0 *           afw   * fw
        aee =  0.125d0 * ( 1.0d0 - afe ) * fe
        ass = -0.125d0 *           afs   * fs
        ann =  0.125d0 * ( 1.0d0 - afn ) * fn

        ap = aw + ae + as + an + aww + aee + ass + ann + df
        !end Quick    

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
        q_art = epsilon1(i,j)* ( pn(i,j)-pn(i,j-1) )/ (y(j)-y(j-1)) + artDivV(i,j)
 
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
!$acc end data
RETURN
END SUBROUTINE RESV1

!#######################################################################################
subroutine upwind_V(um_tau,vm_tau,pn,RV,InvFr2,T,ii,iii,jj,jjj)
    use comum
    ! $acc routine
    implicit none
    integer :: i, j
    integer :: ii, iii, jj, jjj
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: Pn,T
    !REAL(8), DIMENSION(1:imax,1:jmax) :: epsilon1
    real(8) :: InvFr2
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: q_art !artDivU(imax,jmax), artDivV(imax,jmax)
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
    !REAL(8)            :: areau_n(2:imax)            !area n de u
    !REAL(8)            :: areau_s(2:imax)            !area s de u
    !REAL(8)            :: areau_e(1:jmax)            !area e de u
    !REAL(8)            :: areau_w(1:jmax)
    !REAL(8)            :: Re
    !REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    !REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    !REAL(8)            :: b_art         !coeficiente da dissipação artificial
    !REAL(8)            :: liga_poros(imax,jmax)

    !$acc data present(pn,t,liga_poros,rv,areav_s,areav_n,areav_e,areav_w,ym,y,x,xm,epsilon1,vm_tau,um_tau)
    ! $acc data present(re,v_art,residual_v)
    ! $omp parallel do private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dvdydy,dydudx,q_art,v_W,v_E,v_N,v_S,v_P,u_P,ap,aw,aww,ae,aee,an,ann,as,ass)

    ! $acc parallel loop private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dvdydy,dydudx,q_art,v_W,v_E,v_N,v_S,v_P,u_P,ap,aw,aww,ae,aee,an,ann,as,ass)
    !$acc parallel loop private(fw,dw,ds,q_art,v_n,v_e,v_p,u_p,v_w,v_s,an,ae,ap,as,dvdydy,aw,dn,df,de,fs,dydudx,fn,fe)
    do j=jj,jjj
      do i=ii,iii
        fn = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j+1) ) * areav_n(i) / epsilon1(i,j) 
        fs = 0.5d0 * ( vm_tau(i  ,j  ) + vm_tau(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)  
        fe = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)  
        fw = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)  

        df = fe - fw + fn - fs

        Dn = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))

        !CALL upwind(ap,aw,ae,an,as,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
        aw = Dw + MAX(fw , 0.0d0)
        as = Ds + MAX(fs , 0.0d0)
        ae = De + MAX(0.0d0 , -fe)
        an = Dn + MAX(0.0d0 , -fn) 

        ap = aw + ae + as + an + df
 
        v_W = vm_tau(i-1,j  )
        v_E = vm_tau(i+1,j  )
        v_S = vm_tau(i  ,j-1)
        v_N = vm_tau(i  ,j+1)
        v_P = vm_tau(i  ,j  )
        u_P = um_tau(i  ,j  )     !!!!!POROUS NEW

        dvdydy =  areav_n(i) * ( v_N - v_P ) / (ym(j+1)-ym(j  )) &
                 -areav_s(i) * ( v_P - v_S ) / (ym(j  )-ym(j-1)) 
            
        dydudx =  areav_n(i) * (um_tau(i+1,j  )-um_tau(i,j  )) / (xm(i+1)-xm(i)) &
                 -areav_s(i) * (um_tau(i+1,j-1)-um_tau(i,j-1)) / (xm(i+1)-xm(i))
            
        !bulk artificial viscosity term from Ramshaw(1990)
        q_art = epsilon1(i,j)* ( pn(i,j)-pn(i,j-1) )/ (y(j)-y(j-1)) - b_art * (dydudx + dvdydy)
   
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * v_P &
                          + aw * v_W + ae * v_E  &
                          + as * v_S + an * v_N)  &
                          - q_art &
            + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
            -epsilon1(i,j)*  ( v_P/(Re*Darcy_number) + &
            CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
  enddo
enddo
!$acc end data
return
end subroutine upwind_V


!#######################################################################################
subroutine upwind_V1(um_tau,vi,pn,RV,InvFr2,T,ii,iii,jj,jjj)
    use comum
    ! $acc routine
    implicit none
    integer :: i, j
    integer :: ii, iii, jj, jjj
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vi, RV
    REAL(8), DIMENSION(1:imax,1:jmax) :: Pn,T
    !REAL(8), DIMENSION(1:imax,1:jmax) :: epsilon1
    real(8) :: InvFr2
    real(8) :: fw, fe, fs, fn, df, aw, aww, ae, aee, as, ass, an, ann,ap
    REAL(8) :: Dn, Ds, De, Dw
    real(8) :: u_W, u_WW, u_E, u_EE, u_S, u_SS, u_N, u_NN, u_P
    real(8) :: v_W, v_WW, v_E, v_EE, v_S, v_SS, v_N, v_NN, v_P
    real(8) :: q_art !artDivU(imax,jmax), artDivV(imax,jmax)
    real(8) :: dudxdx, dvdydy, dxdvdy, dydudx
    !REAL(8)            :: areau_n(2:imax)            !area n de u
    !REAL(8)            :: areau_s(2:imax)            !area s de u
    !REAL(8)            :: areau_e(1:jmax)            !area e de u
    !REAL(8)            :: areau_w(1:jmax)
    !REAL(8)            :: Re
    !REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    !REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    !REAL(8)            :: b_art         !coeficiente da dissipação artificial
    !REAL(8)            :: liga_poros(imax,jmax)

    !$acc data present(pn,t,liga_poros,rv,areav_s,areav_n,areav_e,areav_w,ym,y,x,xm,epsilon1,vi,um_tau)
    ! $acc data present(re,v_art,residual_v)
    ! $omp parallel do private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dvdydy,dydudx,q_art,v_W,v_E,v_N,v_S,v_P,u_P,ap,aw,aww,ae,aee,an,ann,as,ass)

    ! $acc parallel loop private(i,j,fn,fs,fe,fw,df,Dn,Ds,De,Dw,dvdydy,dydudx,q_art,v_W,v_E,v_N,v_S,v_P,u_P,ap,aw,aww,ae,aee,an,ann,as,ass)
    !$acc parallel loop private(fw,dw,ds,q_art,v_n,v_e,v_p,u_p,v_w,v_s,an,ae,ap,as,dvdydy,aw,dn,df,de,fs,dydudx,fn,fe)
    do j=jj,jjj
      do i=ii,iii
        fn = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j+1) ) * areav_n(i) / epsilon1(i,j) 
        fs = 0.5d0 * ( vi(i  ,j  ) + vi(i  ,j-1) ) * areav_s(i) / epsilon1(i,j)  
        fe = 0.5d0 * ( um_tau(i+1,j  ) + um_tau(i+1,j-1) ) * areav_e(j) / epsilon1(i,j)  
        fw = 0.5d0 * ( um_tau(i  ,j  ) + um_tau(i  ,j-1) ) * areav_w(j) / epsilon1(i,j)  

        df = fe - fw + fn - fs

        Dn = (epsilon1(i,j)/Re) * areav_n(i) / (ym(j+1)-ym(j  ))
        Ds = (epsilon1(i,j)/Re) * areav_s(i) / (ym(j  )-ym(j-1))
        De = (epsilon1(i,j)/Re) * areav_e(j) / (x(i+1)-x(i  ))
        Dw = (epsilon1(i,j)/Re) * areav_w(j) / (x(i  )-x(i-1))

        !CALL upwind(ap,aw,ae,an,as,fw,fe,fn,fs,df,Dw,De,Dn,Ds)
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
            
        dydudx =  areav_n(i) * (um_tau(i+1,j  )-um_tau(i,j  )) / (xm(i+1)-xm(i)) &
                 -areav_s(i) * (um_tau(i+1,j-1)-um_tau(i,j-1)) / (xm(i+1)-xm(i))
            
        !bulk artificial viscosity term from Ramshaw(1990)
        q_art = epsilon1(i,j)* ( pn(i,j)-pn(i,j-1) )/ (y(j)-y(j-1)) - b_art * (dydudx + dvdydy)
   
        RV(i,j) = 1.d0 / (x(i)-x(i-1)) / (y(j)-y(j-1)) * ( - ap * v_P &
                          + aw * v_W + ae * v_E  &
                          + as * v_S + an * v_N)  &
                          - q_art &
            + InvFr2 * ( 1.d0 - 1.d0 / ( (T(i,j)+T(i,j-1)) * 0.5d0 ) ) &
            -epsilon1(i,j)*  ( v_P/(Re*Darcy_number) + &
            CF/((epsilon1(i,j)*Darcy_number)**0.5d0)* v_P*((u_p**2.d0 +v_P**2.d0)**0.5d0))*liga_poros(i,j)
  enddo
enddo
!$acc end data
return
end subroutine upwind_V1


!#######################################################################################
!subroutine solve_V(um,vm,vm_n,um_tau,vm_tau,vm_n_tau,p,T,InvFr2,residual_v,res_v,Rv,vi,areav_e,areav_w,areav_n,areav_s,ym,xm,Re,x,y,liga_poros,epsilon1,artDivV)
subroutine solve_V(um,vm,vm_n,um_tau,vm_tau,vm_n_tau,pn,T,InvFr2,residual_v)
    use comum
    ! $acc routine
    ! $acc routine(RESV)
    ! $acc routine(RESV1)
    ! $acc routine(upwind_v)
    ! $acc routine(upwind_v1)
    ! $acc routine(bcUV)
    implicit none
    integer :: i, j   
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm, vm_n!, vi, RV
    !REAL(8), DIMENSION(2:imax-1,3:jmax-1) :: res_v
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_tau, vm_n_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: pn,T
    !REAL(8), DIMENSION(1:imax  ,1:jmax) :: epsilon1
    REAL(8) :: residual_v
    REAL(8) :: InvFr2
    !REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    !REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    !REAL(8)            :: Re
    !REAL(8)            :: liga_poros(imax,jmax)
    !REAL(8)            :: areav_n(1:imax)            !area n de v
    !REAL(8)            :: areav_s(1:imax)            !area s de v
    !REAL(8)            :: areav_e(2:jmax)            !area e de v
    !REAL(8)            :: areav_w(2:jmax)            !area w de v
    !REAL(8)            :: artDivV(imax,jmax)
    
    !$acc data present(rv,t,invfr2,res_v,liga_poros,areav_e,areav_n,areav_w,xm,y,ym,vm_n_tau,vm_tau,x,areav_s,artdivv,epsilon1,um_tau,vi,pn,vm,um)

    !RV = 0.d0
    !res_v = 0.d0

    CALL RESV(um_tau,vm_tau,pn,T,InvFr2,RV)
    ! $omp parallel do private (i)
    !DO i=2,imax-1
      CALL upwind_V(um_tau,vm_tau,pn,RV,InvFr2,T,2,imax-1,2,2)
    !ENDDO

    ! $omp parallel do private (i)
    !DO i=2,imax-1
      CALL upwind_V(um_tau,vm_tau,pn,RV,InvFr2,T,2,imax-1,jmax,jmax)
    !ENDDO

    ! $omp parallel do private (j)
    !DO j=2,jmax
      CALL upwind_V(um_tau,vm_tau,pn,RV,InvFr2,T,2,2,2,jmax)
    !ENDDO

    ! $omp parallel do private (j)
    !DO j=2,jmax
      CALL upwind_V(um_tau,vm_tau,pn,RV,InvFr2,T,imax-1,imax-1,2,jmax)
    !ENDDO

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j) 
    DO j=3,jmax-1
        DO i=2,imax-1
            res_v(i,j) =( (vm(i,j)-vm_tau(i,j)) +  RV(i,j)*dt) * dtau 
            vi(i,j) = ( vm_tau(i,j) + res_v(i,j) ) 
        ENDDO
    ENDDO

    call bcUV3(um_tau,vi)

    CALL RESV1(um_tau,vi,pn,T,InvFr2,RV)
    ! $omp parallel do private (i)
    !DO i=2,imax-1
      CALL upwind_V1(um_tau,vi,pn,RV,InvFr2,T,2,imax-1,2,2)
    !ENDDO

    ! $omp parallel do private (i)
    !DO i=2,imax-1
      CALL upwind_V1(um_tau,vi,pn,RV,InvFr2,T,2,imax-1,jmax,jmax)
    !ENDDO

    ! $omp parallel do private (j)
    !DO j=2,jmax
      CALL upwind_V1(um_tau,vi,pn,RV,InvFr2,T,2,2,2,jmax)
    !ENDDO

    ! $omp parallel do private (j)
    !DO j=2,jmax
      CALL upwind_V1(um_tau,vi,pn,RV,InvFr2,T,imax-1,imax-1,2,jmax)
    !ENDDO

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j) 
    DO j=3,jmax-1
        DO i=2,imax-1
            res_v(i,j) =( (vm(i,j)-vm_tau(i,j)) +  RV(i,j)*dt) * dtau 
            vi(i,j) =( 0.75d0 * vm_tau(i,j) + 0.25d0 * &
                     ( vi(i,j) + res_v(i,j)) )
        ENDDO
    ENDDO

    call bcUV3(um_tau,vi)

    CALL RESV1(um_tau,vi,pn,T,InvFr2,RV)
    ! $omp parallel do private (i)
    !DO i=2,imax-1
      CALL upwind_V1(um_tau,vi,pn,RV,InvFr2,T,2,imax-1,2,2)
    !ENDDO

    ! $omp parallel do private (i)
    !DO i=2,imax-1
      CALL upwind_V1(um_tau,vi,pn,RV,InvFr2,T,2,imax-1,jmax,jmax)
    !ENDDO

    ! $omp parallel do private (j)
    !DO j=2,jmax
      CALL upwind_V1(um_tau,vi,pn,RV,InvFr2,T,2,2,2,jmax)
    !ENDDO

    ! $omp parallel do private (j)
    !DO j=2,jmax
      CALL upwind_V1(um_tau,vi,pn,RV,InvFr2,T,imax-1,imax-1,2,jmax)
    !ENDDO

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j) 
    DO j=3,jmax-1
        DO i=2,imax-1
            res_v(i,j) =( (vm(i,j)-vm_tau(i,j)) +  RV(i,j)*dt) * dtau 
            vm_n_tau(i,j) = 1.0d0 / 3.0d0 * vm_tau(i,j) + 2.0d0 / 3.0d0 * &
                     ( vi(i,j) + res_v(i,j))
        ENDDO
    ENDDO

    call bcUV4(um_tau,vm_n_tau)
    
    !$acc parallel
    residual_v =  MAXVAL(ABS(res_v)) 
    !$acc end parallel

    !open (550,file='data/resu.dat')    
    !do i = 2,imax-1
    !  do j = 3,jmax-1
    !    write (550,*) X(i) , Y(j) , abs(res_v(i,j))
    !  enddo
    !enddo
    !close(550)
!$acc end data
RETURN
END SUBROUTINE solve_V

!#######################################################################################
!!! Solve Continuity Equation
!subroutine solve_P(p,um_n_tau,vm_n_tau,pn,residual_p,RP,Pi,res_p,areau_e, areau_w, areav_n, areav_s)
SUBROUTINE solve_P(p,um_n_tau,vm_n_tau,pn,residual_p)
    USE comum
    ! $acc routine
    !$acc routine(bcP)
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n_tau
    REAL(8), DIMENSION(1:imax,1:jmax) :: P, Pn!, RP, Pi, res_p
    REAL(8) :: residual_p!, c2, c, max_vel
    REAL(8) :: dudx, dvdy
    !REAL(8)            :: areau_e(1:jmax)            !area e de u
    !REAL(8)            :: areau_w(1:jmax)            !area w de u
    !REAL(8)            :: areav_n(1:imax)            !area n de v
    !REAL(8)            :: areav_s(1:imax)            !area s de v
   
    !$acc data present(rp,vm_n_tau,um_n_tau,residual_p,areau_e,areau_w,areav_n,p,c2,areav_s,pn,pi,res_p)

    !res_p = 0.d0
    ! $acc parallel
    !max_vel = MAXVAL(um_n**2.d0) + MAXVAL(vm_n**2.d0)
    ! $acc end parallel

    ! from AN ARTIFICIAL COMPRESSIBILITY METHOD FOR INCOMPRESSIBLE FLOWS - M. M. Rahman, T. Siikonen
    !c = beta !* SQRT(MAX(max_vel,(0.5d0 * v_i**2.d0))) 
    !c2 = c*c
    c2 = beta !* MAXVAL(vm_n**2.d0)

!!! RALSTON'S METHOD (Second Order Runge-Kutta)
!mass conservation

    ! CALL RESP(um_n,vm_n,p,RP)
    ! $omp parallel do private(i,j,dudx,dvdy)
    !$acc parallel loop collapse(2) private(i,j,dudx,dvdy)
    do j=2,jmax-1
      do i=2,imax-1
        dudx = um_n_tau(i+1,j) * areau_e(j) - um_n_tau(i,j) * areau_w(j)
        dvdy = vm_n_tau(i,j+1) * areav_n(i) - vm_n_tau(i,j) * areav_s(i)
        !AQUII  epsilon
        RP(i,j) = - ( dudx + dvdy )    
      enddo
    enddo

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j)
    DO j=2,jmax-1
      DO i=2,imax-1      
        pi(i,j) = p(i,j) + dtau * RP(i,j)* c2 
      ENDDO
    ENDDO

    CALL bcP1(pi)

    !CALL RESP(um_n,vm_n,pi,RP)
    ! $omp parallel do private(i,j,dudx,dvdy)
    !$acc parallel loop collapse(2) private(i,j,dudx,dvdy)
    do j=2,jmax-1
      do i=2,imax-1
        dudx = um_n_tau(i+1,j) * areau_e(j) - um_n_tau(i,j) * areau_w(j)
        dvdy = vm_n_tau(i,j+1) * areav_n(i) - vm_n_tau(i,j) * areav_s(i)
        RP(i,j) = - ( dudx + dvdy )    
      enddo
    enddo

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j)
    DO j=2,jmax-1
      DO i=2,imax-1          
        pi(i,j) = 0.75d0 * p(i,j) + 0.25d0 * (pi(i,j) + dtau * RP(i,j)* c2)
      ENDDO
    ENDDO

    CALL bcP1(pi)

    !CALL RESP(um_n,vm_n,pi,RP)
    ! $omp parallel do private(i,j,dudx,dvdy)
    !$acc parallel loop collapse(2) private(i,j,dudx,dvdy)
    do j=2,jmax-1
      do i=2,imax-1
        dudx = um_n_tau(i+1,j) * areau_e(j) - um_n_tau(i,j) * areau_w(j)
        dvdy = vm_n_tau(i,j+1) * areav_n(i) - vm_n_tau(i,j) * areav_s(i)
        RP(i,j) = - ( dudx + dvdy )    
      enddo
    enddo
    
    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j)
    DO j=2,jmax-1
      DO i=2,imax-1
        res_p(i,j) = dtau * RP(i,j) * c2
        pn(i,j) = 1.0d0 / 3.0d0 * p(i,j) + 2.0d0 / 3.0d0 * (pi(i,j) + res_p(i,j))
      ENDDO
    ENDDO

    CALL bcP(pn)

    !$acc parallel
    residual_p = MAXVAL(ABS(res_p))
    !$acc end parallel
!$acc end data
RETURN
END SUBROUTINE solve_P

!#######################################################################################
!!! Solve Mixture Fraction
!subroutine solve_Z(um_n_tau,vm_n_tau,T,T_n_tau,T_tau,res_Z,RZ,zi,areau_e,areau_w,areav_n,areav_s,ym,xm,Pe,x,y,liga_poros,dtau)
SUBROUTINE solve_Z(um_n_tau,vm_n_tau,T,T_n_tau,T_tau)
    USE comum
    ! $acc routine
    ! $acc routine(resZ)
    !$acc routine(bcZ)
    IMPLICIT NONE
    INTEGER :: i, j
    REAL(8), DIMENSION(1:imax+1,1:jmax  ) :: um_n_tau
    REAL(8), DIMENSION(1:imax  ,1:jmax+1) :: vm_n_tau
    !REAL(8), DIMENSION(1:imax,1:jmax) :: RZ, Zi
    REAL(8), DIMENSION(1:imax,1:jmax) :: T_tau, T, T_n_tau
    !REAL(8), DIMENSION(2:imax-1,2:jmax-1) :: res_Z
    REAL(8) :: dZudx, dZvdy

    REAL(8) :: Dw,De,Ds,Dn, Dp
    REAL(8) :: Zw,Ze,Zs,Zn, Zp
    !REAL(8)            :: x(1:imax), y(1:jmax)       !malha principal
    !REAL(8)            :: xm(1:imax+1), ym(1:jmax+1) !malha deslocada
    !REAL(8)            :: areau_e(1:jmax)            !area e de u
    !REAL(8)            :: areau_w(1:jmax)            !area w de u
    !REAL(8)            :: areav_n(1:imax)            !area n de v
    !REAL(8)            :: areav_s(1:imax)            !area s de v
    !REAL(8)            :: Pe,dtau
    !REAL(8)            :: liga_poros(imax,jmax)

    !$acc data present(areav_s,areav_n,x,vm_n_tau,um_n_tau,xm,areau_w,areau_e,y,t_tau,t_n_tau,rz,res_z,liga_poros,zi,ym,t)

    !REAL(8), DIMENSION(1:imax,1:jmax ) :: epsilon1
    !RALSTON'S METHOD (Second Order Runge-Kutta)
    !CALL RESZ(um_n,vm_n,Z_tau,RZ)
    ! $omp parallel do private(i,j,dZudx,dZvdy,De,Dw,Dn,Ds,Dp) 
 
    !$acc parallel loop collapse(2) private(i,j,dZudx,dZvdy,De,Dw,Dn,Ds,Dp) 
    do j=2,jmax-1
      do i=2,imax-1
        dZudx = 0.5d0 * (T(i+1,j)+T(i,j)) * um_n_tau(i+1,j) * areau_e(j) &
              - 0.5d0 * (T(i-1,j)+T(i,j)) * um_n_tau(i  ,j) * areau_w(j)

        dZvdy = 0.5d0 * (T(i,j+1)+T(i,j)) * vm_n_tau(i,j+1) * areav_n(i) &
              - 0.5d0 * (T(i,j-1)+T(i,j)) * vm_n_tau(i,j  ) * areav_s(i)

        De = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i+1)-x(i  ))  
        Dw = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i  )-x(i-1))  
        Dn = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j+1)-y(j  ))  
        Ds = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j  )-y(j-1))  
        Dp = De + Dw + Dn + Ds 

        RZ(i,j) = 1.d0 / (xm(i+1)-xm(i)) / (ym(j+1)-ym(j)) *&
            (-Dp*T(i,j  ) + De*T(i+1,j) + Dw*T(i-1,j) + Dn*T(i,j+1) + Ds*T(i,j-1) - (1.d0-liga_poros(i,j))* (dZudx + dZvdy) )!/epsilon1(i,j))
      enddo
    enddo

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j) 
    DO j=2,jmax-1
      DO i=2,imax-1
        res_Z(i,j) =( (T(i,j)-T_tau(i,j)) +  RZ(i,j)*dt) * dtau 
        Zi(i,j) = T_tau(i,j) + res_Z(i,j)
      ENDDO
    ENDDO

    CALL bcZ2(Zi)

    !CALL RESZ(um_n,vm_n,Zi,RZ)
    ! $omp parallel do private(i,j,dZudx,dZvdy,De,Dw,Dn,Ds,Dp) 
    !$acc parallel loop collapse(2) private(i,j,dZudx,dZvdy,De,Dw,Dn,Ds,Dp) 
    do j=2,jmax-1
      do i=2,imax-1
        dZudx = 0.5d0 * (T(i+1,j)+T(i,j)) * um_n_tau(i+1,j) * areau_e(j) &
              - 0.5d0 * (T(i-1,j)+T(i,j)) * um_n_tau(i  ,j) * areau_w(j)

        dZvdy = 0.5d0 * (T(i,j+1)+T(i,j)) * vm_n_tau(i,j+1) * areav_n(i) &
              - 0.5d0 * (T(i,j-1)+T(i,j)) * vm_n_tau(i,j  ) * areav_s(i)

        De = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i+1)-x(i  ))
        Dw = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i  )-x(i-1))
        Dn = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j+1)-y(j  ))
        Ds = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j  )-y(j-1))
        Dp = De + Dw + Dn + Ds 

        RZ(i,j) = 1.d0 / (xm(i+1)-xm(i)) / (ym(j+1)-ym(j)) *&
            (-Dp*T(i,j  ) + De*T(i+1,j) + Dw*T(i-1,j) + Dn*T(i,j+1) + Ds*T(i,j-1) - (1.d0-liga_poros(i,j))* (dZudx + dZvdy) )!/epsilon1(i,j))
      enddo
    enddo

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j) 
    DO j=2,jmax-1
      DO i=2,imax-1
        res_Z(i,j) =( (T(i,j)-T_tau(i,j)) +  RZ(i,j)*dt) * dtau 
        Zi(i,j) = 0.75d0 * T_tau(i,j) + 0.25d0 * (Zi(i,j) + res_Z(i,j))
      ENDDO
    ENDDO

    CALL bcZ2(Zi)

    !CALL RESZ(um_n,vm_n,Zi,RZ)
    ! $omp parallel do private(i,j, dZudx, dZvdy, De, Dw, Dn, Ds, Ze,Zw,Zn,Zs,Zp, Dp) 
    !$acc parallel loop collapse(2) private(i,j,dZudx,dZvdy,De,Dw,Dn,Ds,Dp) 
    do j=2,jmax-1
      do i=2,imax-1
        dZudx = 0.5d0 * (T(i+1,j)+T(i,j)) * um_n_tau(i+1,j) * areau_e(j) &
              - 0.5d0 * (T(i-1,j)+T(i,j)) * um_n_tau(i  ,j) * areau_w(j)

        dZvdy = 0.5d0 * (T(i,j+1)+T(i,j)) * vm_n_tau(i,j+1) * areav_n(i) &
              - 0.5d0 * (T(i,j-1)+T(i,j)) * vm_n_tau(i,j  ) * areav_s(i)

        De = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i+1)-x(i  ))  
        Dw = (ym(j+1)-ym(j)) * (1.d0/Pe) / (x(i  )-x(i-1))  
        Dn = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j+1)-y(j  ))  
        Ds = (xm(i+1)-xm(i)) * (1.d0/Pe) / (y(j  )-y(j-1))  
        Dp = De + Dw + Dn + Ds 

        RZ(i,j) = 1.d0 / (xm(i+1)-xm(i)) / (ym(j+1)-ym(j)) *&
            (-Dp*T(i,j  ) + De*T(i+1,j) + Dw*T(i-1,j) + Dn*T(i,j+1) + Ds*T(i,j-1) - (1.d0-liga_poros(i,j))* (dZudx + dZvdy) )!/epsilon1(i,j))
      enddo
    enddo

    ! $omp parallel do private(i,j) 
    !$acc parallel loop collapse(2) private(i,j) 
    DO j=2,jmax-1
      DO i=2,imax-1
        res_Z(i,j) =( (T(i,j)-T_tau(i,j)) +  RZ(i,j)*dt) * dtau 
        T_n_tau(i,j) = 1.0d0 / 3.0d0 * T_tau(i,j) + 2.0d0 / 3.0d0 * (Zi(i,j) + res_Z(i,j))
        ENDDO
    ENDDO

    CALL bcZ1(T_n_tau)
!$acc end data
RETURN
END SUBROUTINE solve_Z

