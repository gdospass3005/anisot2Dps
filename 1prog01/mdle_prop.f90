MODULE mdle_prop

CONTAINS
  SUBROUTINE prop_anisot(Vx,Vz,Sxx,Sxz,Szz,&
       C11,C13,C33,C44,rhox,dx,dt,Nx,Nz,&
       trig_x,ifax_x,nfac_x,rkx,trig_z,ifax_z,nfac_z,rkz)
    USE mdle_fft
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: Vx,Vz,Sxx,Sxz,Szz
    REAL, DIMENSION(:,:), INTENT(IN)    :: C11,C13,C33,C44
    REAL, DIMENSION(:,:), INTENT(IN)    :: rhox
    REAL,    INTENT(IN) :: dx,dt
    INTEGER, INTENT(IN) :: Nx,Nz
    REAL, DIMENSION(2*Nx), INTENT(IN)  :: trig_x
    REAL, DIMENSION(2*Nz), INTENT(IN)  :: trig_z
    INTEGER, DIMENSION(Nx), INTENT(IN)  :: ifax_x
    INTEGER, DIMENSION(Nz), INTENT(IN)  :: ifax_z
    INTEGER, INTENT(IN)  :: nfac_x,nfac_z
    REAL, DIMENSION(Nx), INTENT(IN)  :: rkx
    REAL, DIMENSION(Nz), INTENT(IN)  :: rkz

    REAL, DIMENSION(Nx) :: aux1
    REAL, DIMENSION(Nz) :: aux2
    INTEGER :: i,j
    REAL :: l2dx

    l2dx=1./(2.*dx);

    ! Space derivatives of stresses => velocities
    DO j=1,Nz
       aux1=Sxx(:,j)
       call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
       Vx(:,j) = Vx(:,j) + (dt * rhox(:,j) * aux1)
    END DO

    DO i=1,Nx
       aux2=Sxz(i,:)
       call ddx(aux2,Nz,trig_z,ifax_z,nfac_z,rkz)
       Vx(i,:) = Vx(i,:) + (dt * rhox(i,:) * aux2)
    END DO

    DO j=1,Nz
       aux1=Sxz(:,j)
       call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
       Vz(:,j) = Vz(:,j) + (dt * rhox(:,j) * aux1)
    END DO

    DO i=1,Nx
       aux2=Szz(i,:)
       call ddx(aux2,Nz,trig_z,ifax_z,nfac_z,rkz)
       Vz(i,:) = Vz(i,:) + (dt * rhox(i,:) * aux2)
    END DO

    ! Space derivatives of velocities => stresses
    DO j=1,Nz
       aux1=Vx(:,j)
       call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
       Sxx(:,j) = Sxx(:,j) + (dt * C11(:,j) * aux1)
       Szz(:,j) = Szz(:,j) + (dt * C13(:,j) * aux1)
    END DO

    DO i=1,Nx
       aux2=Vz(i,:)
       call ddx(aux2,Nz,trig_z,ifax_z,nfac_z,rkz)
       Sxx(i,:) = Sxx(i,:) + (dt * C13(i,:) * aux2)
       Szz(i,:) = Szz(i,:) + (dt * C33(i,:) * aux2)
    END DO

    DO j=1,Nz
       aux1=Vz(:,j)
       call ddx(aux1,Nx,trig_x,ifax_x,nfac_x,rkx)
       Sxz(:,j) = Sxz(:,j) + (dt * C44(:,j) * aux1)
    END DO

    DO i=1,Nx
       aux2=Vx(i,:)
       call ddx(aux2,Nz,trig_z,ifax_z,nfac_z,rkz)
       Sxz(i,:) = Sxz(i,:) + (dt * C44(i,:) * aux2)
    END DO

  END SUBROUTINE prop_anisot

END MODULE mdle_prop



