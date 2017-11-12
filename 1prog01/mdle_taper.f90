MODULE mdle_taper

CONTAINS

  SUBROUTINE bt_exp_create(taper,nb,F)
    IMPLICIT NONE
    REAL :: F
    INTEGER :: nb
    REAL, DIMENSION(nb) :: taper
    INTEGER :: i

    DO i=1,nb
       taper(i) = exp( -(F*(REAL(nb) - REAL(i) ))**2 )
    END DO

  END SUBROUTINE bt_exp_create


  SUBROUTINE bt_apply_multiple(Vx,Vz,Sxx,Sxz,Szz,&
       Nx,Nz,nb,taper,fsf)
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: Vx,Vz,Sxx,Sxz,Szz
    INTEGER, INTENT(IN) :: Nx,Nz,nb,fsf
    REAL, DIMENSION(nb), INTENT(IN) :: taper

    CALL bt_apply(Vx,Nx,Nz,nb,taper,fsf)
    CALL bt_apply(Vz,Nx,Nz,nb,taper,fsf)
    CALL bt_apply(Sxx,Nx,Nz,nb,taper,fsf)
    CALL bt_apply(Sxz,Nx,Nz,nb,taper,fsf)
    CALL bt_apply(Szz,Nx,Nz,nb,taper,fsf)

  END SUBROUTINE bt_apply_multiple


  SUBROUTINE bt_apply(pp,Nx,Nz,nb,taper,fsf)
    IMPLICIT NONE
    INTEGER :: NX,NZ,nb,fsf
    REAL, DIMENSION(nb) :: taper
    REAL, DIMENSION(NX,NZ) :: pp
    INTEGER :: i

    DO i=1,NZ
       pp(1:nb,i) = pp(1:nb,i) * taper
       pp(NX:NX-nb+1:-1,i) =  pp(NX:NX-nb+1:-1,i) * taper
    END DO

    DO i=1,NX
       IF (fsf.ne.1) THEN
          pp(i,1:nb:1) = pp(i,1:nb:1) * taper
       END IF
       pp(i,NZ:NZ-nb+1:-1) =  pp(i,NZ:NZ-nb+1:-1) * taper
    END DO

  END SUBROUTINE bt_apply

END MODULE mdle_taper



