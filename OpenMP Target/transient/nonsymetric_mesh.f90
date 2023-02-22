SUBROUTINE mesh
    USE comum
    IMPLICIT NONE
    INTEGER :: i, j, ii, jj
    !parametros de refinamento, max = 2
    REAL(8) :: etax, etay, sx, sy


    DO i=1,imax                              !o termo abaixo tira a simetria da malha
        X(i) = ( (dfloat(i) - 1.d0) * dx_c )  - Lhori/2.d0
    ENDDO

    DO j=1,jmax
        Y(j) = ( (dfloat(j) - 1.d0) * dx_c ) - y_down
    ENDDO

!############## APLICA O REFINAMENTO NA MALHA CONSTANTE #######################
    !aplica a função de stretching para x
!    DO i = 1,imax
!            etax = ( X(i) - X(imax) ) / ( X(1) - X(imax) )
!            sx = Px_grid * etax + (1.d0 - Px_grid) &
!                 * (1.d0 - ( (tanh(Q_grid * (1.d0 - etax)) ) / tanh(Q_grid) ) )
!            X(i) = X(imax) - sx * (X(imax) - X(1) )
!    ENDDO

    DO i = imax/2+1,imax
            etax = ( X(i) - X(imax) ) / ( X(imax/2+1) - X(imax) )
            sx = Px_grid * etax + (1.d0 - Px_grid) &
                 * (1.d0 - ( (tanh(Q_grid * (1.d0 - etax)) ) / tanh(Q_grid) ) )
            X(i) = X(imax) - sx * (X(imax) - X(imax/2+1) )
    ENDDO

    DO i = imax/2+1,1,-1
            etax = ( X(i) + X(imax) ) / ( X(imax/2+1) + X(imax) )
            sx = Px_grid * etax + (1.d0 - Px_grid) &
                 * (1.d0 - ( (tanh(Q_grid * (1.d0 - etax)) ) / tanh(Q_grid) ) )
            X(i) = -X(imax) - sx * (-X(imax) - X(imax/2+1) )
    ENDDO

    !aplica a função de stretching para y acima do cilindro
    DO j = int(y_down/dx_c)+1,jmax
        etay = ( Y(j) - Y(jmax) ) / ( Y(int(y_down/dx_c)+1) - Y(jmax) )
        sy = Py_grid * etay + (1.d0 - Py_grid) &
             * (1.d0 - ( (tanh(Q_grid * (1.d0 - etay)) ) / tanh(Q_grid) ) )
        Y(j) = Y(jmax) - sy * (Y(jmax) - Y(int(y_down/dx_c)+1) )
    ENDDO
    

    DO j = int(y_down/dx_c)+1,1,-1
        etay = ( Y(j) + Y(jmax) ) / ( Y(int(y_down/dx_c)+1) + Y(jmax) )
        sy = Py_grid * etay + (1.d0 - Py_grid) &
             * (1.d0 - ( (tanh(Q_grid * (1.d0 - etay)) ) / tanh(Q_grid) ) )
        Y(j) = -Y(jmax) - sy * (-Y(jmax) - Y(int(y_down/dx_c)+1) )
    ENDDO


!##############CALCULA OS PONTOS MEDIOS DA MALHA XM E YM ######################

    xm(1) = -(x(2)+x(1))*0.5d0
    xm(imax+1) = (x(imax)+x(imax-1))*0.5d0 + (x(imax)-x(imax-1))

    DO i= 2,imax
        xm(i) = (x(i)+x(i-1))*0.5d0
    ENDDO

    ym(1) = y(1) - (y(2)-y(1))*0.5d0
    ym(jmax+1) = (y(jmax)+y(jmax-1))*0.5d0 + (y(jmax)-y(jmax-1))

    DO j= 2,jmax
        ym(j) = (y(j)+y(j-1))*0.5d0 
    ENDDO

!###############CALCULA OS VOLUMES DE U, V E P ###################################

   DO i=2,imax
       DO j=1,jmax
            vol_u(i,j) = (x(i)-x(i-1)) * (ym(j+1)-ym(j)) 
       ENDDO
   ENDDO

   DO i=1,imax
       DO j=2,jmax
            vol_v(i,j) = (y(j)-y(j-1)) * (xm(i+1)-xm(i)) 
       ENDDO
   ENDDO


   DO i=1,imax
       DO j=2,jmax
            vol_p(i,j) = (ym(j+1)-ym(j)) * (xm(i+1)-xm(i)) 
       ENDDO
   ENDDO
!###############CALCULA AS AREAS DAS FACES W,E,N,S DE U E V ###################

    DO i=2,imax
        areau_n(i) = x(i) - x(i-1)
        areau_s(i) = x(i) - x(i-1)
    ENDDO

    DO j=1,jmax
        areau_e(j) = ym(j+1) - ym(j)
        areau_w(j) = ym(j+1) - ym(j)
    ENDDO


    DO i=1,imax
        areav_n(i) = xm(i+1) - xm(i)
        areav_s(i) = xm(i+1) - xm(i)
    ENDDO

    DO j=2,jmax
        areav_e(j) = y(j) - y(j-1)
        areav_w(j) = y(j) - y(j-1)
    ENDDO

!##########################################################

    DO j=2,jmax
        dy(j) = y(j)-y(j-1)
    ENDDO

    DO i=2,imax
        dx(i) = x(i)-x(i-1)
    ENDDO

            WRITE(*,*) 'dx =',maxval(dx) , '[mm]'
            WRITE(*,*) 'dy =',maxval(dy) , '[mm]'

    OPEN (550,file='data/mesh.dat')

        WRITE (550,*) imax, jmax , Hvert

    CLOSE(550)

    OPEN (1,file="data/grid.dat")

    DO j=1,jmax,2 !plotar na direção de i
        DO i=1,imax
            WRITE(1,*) x(i),y(j)
        END DO

        IF (j<jmax) THEN
            jj=j+1
        DO i=imax,1,-1
            WRITE(1,*) x(i),y(jj)
        END DO
        END IF
    END DO

    DO i=1,imax,2 !plotar na direção de j
        DO j=jmax,1,-1
            WRITE(1,*) x(i),y(j)
        END DO

        IF (i<imax) THEN
            ii=i+1
            DO j=1,jmax
            WRITE(1,*) x(ii),y(j)
            END DO
        END IF
    END DO
    CLOSE(1)

 
    OPEN (110,file="data/grid_u.dat")
        DO i=1,imax+1
            DO j=1,jmax
                WRITE(110,*) xm(i), y(j)
            ENDDO
        ENDDO

    CLOSE(110)

    OPEN (110,file="data/grid_v.dat")
        DO i=1,imax
            DO j=1,jmax+1
                WRITE(110,*) x(i), ym(j)
            ENDDO
        ENDDO

    CLOSE(110)

!------------- cilindro r=1 ------------------------------------
    DO i = 1,imax
        DO j = 1, jmax
        IF (x(i) * x(i) + y(j) * y(j) <= ao/L_c * ao/L_c ) THEN 

            flag(i,j) = C_I

        ELSE 
           flag(i,j) = C_F
        END IF

        END DO
    END DO

!------------bloco L=1 ----------------------------------------
!    DO i = 1,imax
!        DO j = 1, jmax
!        IF (x(i) .LE. 0.5d0 .AND. y(j) .LE. 0.5d0 .AND. y(j) .GE. -0.5d0 ) THEN 
!
!            flag(i,j) = C_I
!
!        ELSE 
!           flag(i,j) = C_F
!        END IF
!
!        END DO
!    END DO


!--------- sem objeto no dominio -------------
!    DO i = 1,imax
!        DO j = 1, jmax
!           flag(i,j) = C_F
!        END DO
!    END DO
!

    DO i = 2, imax-1
        DO j = 2, jmax-1
            IF (flag(i,j) .EQ. C_I .AND. flag(i-1,j  ) .EQ. C_F &
            .or.flag(i,j) .EQ. C_I .AND. flag(i  ,j-1) .EQ. C_F &
            .or.flag(i,j) .EQ. C_I .AND. flag(i+1,j  ) .EQ. C_F&
            .or.flag(i,j) .EQ. C_I .AND. flag(i  ,j+1) .EQ. C_F) THEN
        
            flag(i,j) = C_B

            ENDIF
        ENDDO
    ENDDO

    i=1
    DO j=2, jmax-1
        IF (flag(i,j).EQ.C_I.AND.flag(i,j-1).EQ.C_F &
        .or.flag(i,j).EQ.C_I.AND.flag(i,j+1).EQ.C_F) THEN

        flag(i,j) = C_B

        ENDIF
    ENDDO


    DO i = 1, imax-1
        DO j = 1, jmax-1
            IF (flag(i,j).EQ.C_F.AND.flag(i-1,j  ).EQ.C_B &
            .or.flag(i,j).EQ.C_F.AND.flag(i  ,j-1).EQ.C_B &
            .or.flag(i,j).EQ.C_F.AND.flag(i+1,j  ).EQ.C_B &
            .or.flag(i,j).EQ.C_F.AND.flag(i  ,j+1).EQ.C_B) THEN

            flag(i,j) = C_BS

            ENDIF
        ENDDO
    ENDDO



!on/off for porosity
    DO i = 1,imax
        DO j = 1, jmax
        IF (flag(i,j) .NE. C_F) THEN
        
          epsilon1(i,j) =   porosidade       
          liga_poros(i,j) = 1.d0
            ELSE 
          epsilon1(i,j) = 1.d0
          liga_poros(i,j) = 0.d0
         END IF
write(*,*) i,j,epsilon1(i,j), liga_poros(i,j)
        END DO
    END DO




 


    OPEN (1,file="data/grid_droplet.dat")
        DO i = 1, imax
            DO j = 1, jmax

            IF (flag(i,j).EQ.C_I) THEN 
              WRITE(1,*) x(i),y(j)
            ENDIF

            ENDDO
        ENDDO
    CLOSE(1)

    OPEN (1,file="data/grid_boundary.dat")
        DO i = 1, imax
            DO j = 1, jmax

            IF (flag(i,j).EQ.C_B) THEN 
              WRITE(1,*) x(i),y(j)
            ENDIF

            ENDDO
        ENDDO
    CLOSE(1)

    OPEN (2,file="data/grid_boundary_side.dat")
        DO i = 1, imax
            DO j = 1, jmax

            IF (flag(i,j).EQ.C_BS) THEN 
              WRITE(2,*) x(i),y(j)
            ENDIF

            ENDDO
        ENDDO
    CLOSE(2)



RETURN
END SUBROUTINE mesh

