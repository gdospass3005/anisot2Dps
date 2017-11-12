PROGRAM anisot

  USE mdle_source
  USE mdle_prop
  USE mdle_taper
  USE mdle_io_utils                                                  
  USE mdle_fft

  IMPLICIT NONE
  REAL, DIMENSION(:,:), ALLOCATABLE :: Vx,Vz,Sxx,Sxz,Szz

  INTEGER :: Nx,Nz
  REAL    :: dx,dz,dt,fpeak
  INTEGER :: itmax
  INTEGER :: lx, lz
  REAL, DIMENSION(:,:), ALLOCATABLE :: C11,C13,C33,C44
  REAL, DIMENSION(:,:), ALLOCATABLE :: rhox

  REAL, DIMENSION(:),   ALLOCATABLE :: trs
  REAL    :: tdelay

  REAL, DIMENSION(:), ALLOCATABLE :: taper  
  REAL :: F
  INTEGER :: nb,fsf

  REAL    :: t
  INTEGER :: it,i,j

  REAL, DIMENSION(:,:), ALLOCATABLE :: csg_Vx,csg_Vz,csg_P,csg_S
  INTEGER :: nphones,npmin_x,npmin_z,dnp_x,dnp_z
  INTEGER :: n1,n2,n3,n4,n5,n6,n7,n8
  INTEGER :: isnap,nsnaps,snapmin,dsnap

  REAL, DIMENSION(:), ALLOCATABLE  :: trig_x
  REAL, DIMENSION(:), ALLOCATABLE  :: trig_z
  INTEGER, DIMENSION(:), ALLOCATABLE  :: ifax_x
  INTEGER, DIMENSION(:), ALLOCATABLE  :: ifax_z
  INTEGER  :: nfac_x,nfac_z
  REAL, DIMENSION(:), ALLOCATABLE  :: rkx
  REAL, DIMENSION(:), ALLOCATABLE  :: rkz

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  PRINT*,'anisot_ps - Modeling seismic waves in anisotropic media'
  PRINT*,'            (Pseudospectral method)'

  ! Input
  CALL inputdata(Nx,Nz,dx,dt,fpeak,itmax,lx,lz,&
       nphones,npmin_x,npmin_z,dnp_x,dnp_z,nsnaps,snapmin,dsnap,&
       nb,F,fsf)

  dz=dx

  ALLOCATE(C11(Nx,Nz),C13(Nx,Nz),C33(Nx,Nz),C44(Nx,Nz))
  ALLOCATE(rhox(Nx,Nz))

  CALL inputmodel(Nx,Nz,C11,C13,C33,C44,rhox)
  C11 = C11*1.e10
  C13 = C13*1.e10
  C33 = C33*1.e10
  C44 = C44*1.e10

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ALLOCATE(trs(itmax))
  CALL source_init(trs,dt,itmax,fpeak,tdelay)

  ALLOCATE(taper(nb))
  CALL bt_exp_create(taper,nb,F)

  ! Initialize fft vectors
  ALLOCATE(trig_x(2*Nx),ifax_x(Nx),rkx(Nx))
  ALLOCATE(trig_z(2*Nz),ifax_z(Nz),rkz(Nz))
  call fft_init(NX,trig_x,ifax_x,nfac_x,dx,rkx)
  call fft_init(NZ,trig_z,ifax_z,nfac_z,dz,rkz)

  ALLOCATE(csg_Vx(itmax,nphones),csg_Vz(itmax,nphones))
  ALLOCATE(csg_P(itmax,nphones),csg_S(itmax,nphones))

  n1=31; n2=32; n3=33; n4=34; n5=35; n6=36; n7=37; n8=38
  OPEN(n1,FILE="csg_Vx.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  OPEN(n2,FILE="csg_Vz.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)

  OPEN(n3,FILE="snapshots_Vx.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)
  OPEN(n4,FILE="snapshots_Vz.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

  OPEN(n5,FILE="csg_P.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
  OPEN(n6,FILE="snapshots_P.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
       ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

   OPEN(n7,FILE="csg_S.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
        ACTION='WRITE', FORM='UNFORMATTED', RECL=itmax*4)
   OPEN(n8,FILE="snapshots_S.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
        ACTION='WRITE', FORM='UNFORMATTED', RECL=Nz*4)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ALLOCATE(Vx(Nx,Nz),Vz(Nx,Nz),Sxx(Nx,Nz),Sxz(Nx,Nz),Szz(Nx,Nz))
  Sxx=0.; Sxz=0.; Szz=0.; Vx=0.; Vz=0.
  rhox=1./rhox;
  isnap=0

  ! Beggining time loop
  DO it=1,itmax
     t=it*dt
     OPEN(77,FILE="status.txt",STATUS='UNKNOWN',ACTION='WRITE')
     WRITE(77,*) &
          'anisot - Modeling seismic waves in anisotropic media'
     WRITE(77,*) &
          'Iteracao',it,'/',itmax,' Vx(lx+5,lz+5)=',Vx(lx+5,lz+5)
     CLOSE(77)

     ! Insert source function 
     Sxx(lx,lz)=Sxx(lx,lz)+trs(it)
     Szz(lx,lz)=Szz(lx,lz)+trs(it)

     ! Do one time step
     CALL prop_anisot(Vx,Vz,Sxx,Sxz,Szz,&
          C11,C13,C33,C44,rhox,dx,dt,Nx,Nz,&
             trig_x,ifax_x,nfac_x,rkx,trig_z,ifax_z,nfac_z,rkz)

     PRINT*,'it',it,'/',itmax,' Vx(lx+5,lz+5)=',Vx(lx+5,lz+5)

     ! Output
     CALL save_shotgather_n_snapshots(csg_Vx,csg_Vz,&
          csg_P,csg_S,Vx,Vz,&
          Sxx,Sxz,Szz,nphones,npmin_x,npmin_z,dnp_x,dnp_z,&
          isnap,nsnaps,snapmin,dsnap,n1,n2,n3,n4,n5,n6,n7,n8,&
          Nx,Nz,it,itmax)

     ! Boundary taper
     CALL bt_apply_multiple(Vx,Vz,Sxx,Sxz,Szz,Nx,Nz,nb,taper,fsf)

  END DO  !End time loop

  PRINT*,'Successful run.'

END PROGRAM anisot
MODULE mdle_fft

INTERFACE
SUBROUTINE ffttrig(N,trig)
real, DIMENSION(2*N), INTENT(OUT) :: trig
INTEGER, INTENT(IN) :: N
END SUBROUTINE ffttrig
SUBROUTINE factork(ifax,nfac,N)
INTEGER, INTENT(IN) :: N
INTEGER, DIMENSION(N), INTENT(OUT) :: ifax
INTEGER, INTENT(OUT) :: nfac
END SUBROUTINE factork
SUBROUTINE fft(data,cccc,N,trig,ifax,nfac,skip,isign)
REAL, DIMENSION(2*N), INTENT(INOUT) :: data
REAL, DIMENSION(2*N), INTENT(INOUT) :: cccc
REAL, DIMENSION(2*N), INTENT(IN) :: trig
INTEGER, INTENT(IN) :: skip
INTEGER, INTENT(IN) :: N,isign
INTEGER, DIMENSION(N), INTENT(IN) :: ifax
INTEGER, INTENT(IN) :: nfac
END SUBROUTINE fft
END INTERFACE

CONTAINS

SUBROUTINE fft_init(N,trig,ifax,nfac,dx,rkx)
INTEGER, INTENT(IN) :: N
REAL, DIMENSION(2*N), INTENT(OUT) :: trig
INTEGER, DIMENSION(N), INTENT(OUT) :: ifax
INTEGER, INTENT(OUT) :: nfac
REAL, INTENT(IN) :: dx
REAL, DIMENSION(N), INTENT(OUT) :: rkx
REAL :: PI=3.141592653589793238462643383279502884197
REAL :: dkx
INTEGER :: inyq_kx

call ffttrig(N,trig)
call factork(ifax,nfac,N)

! Spatial frequency sampling interval
dkx = 2*PI/(N*dx)
! Location of the Nyquist frequency on the frequency vector rkx
inyq_kx = N/2 + 1
! Nyquist spatial frequency
!kx_nyq = inyq_kx * dkx

! Frequency vector
do i=1,inyq_kx
  rkx(i) = (i-1)*dkx
end do
do i=inyq_kx+1,N
  rkx(i) = -rkx(N+2-i)
end do

END SUBROUTINE fft_init

SUBROUTINE ddx(data,N,trig,ifax,nfac,rkx)
IMPLICIT NONE
REAL, DIMENSION(N), INTENT(INOUT) :: data
INTEGER, INTENT(IN) :: N
REAL, DIMENSION(2*N), INTENT(IN) :: trig
INTEGER, DIMENSION(N), INTENT(IN) :: ifax
INTEGER, INTENT(IN) :: nfac
REAL, DIMENSION(N), INTENT(IN) :: rkx
REAL, DIMENSION(2*N) :: datadouble
REAL, DIMENSION(2*N) :: cccc
COMPLEX, DIMENSION(N) :: cdata
COMPLEX :: j

datadouble(1:2*N:2) = data
datadouble(2:2*N:2) = 0.0
call fft(datadouble,cccc,N,trig,ifax,nfac,1,1)
cdata=cmplx(datadouble(1:2*N:2),datadouble(2:2*N:2))
j = cmplx(0.,1.)
cdata=rkx*j*cdata
datadouble(1:2*N:2)=real(cdata)
datadouble(2:2*N:2)=aimag(cdata)
call fft(datadouble,cccc,N,trig,ifax,nfac,1,-1)

data = datadouble(1:2*N:2)/N

END SUBROUTINE ddx

END MODULE mdle_fft
MODULE mdle_io_utils

CONTAINS

  SUBROUTINE inputdata(Nx,Nz,dx,dt,fpeak,itmax,lx,lz,&
       nphones,npmin_x,npmin_z,dnp_x,dnp_z,nsnaps,snapmin,dsnap,&
       nb,F,fsf)
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: Nx,Nz,itmax,lx,lz
    REAL,    INTENT(OUT) :: dx,dt,fpeak
    INTEGER, INTENT(OUT) :: nphones,npmin_x,npmin_z,dnp_x,dnp_z
    INTEGER, INTENT(OUT) :: nsnaps,snapmin,dsnap
    INTEGER, INTENT(OUT) :: nb,fsf
    REAL,    INTENT(OUT) :: F
 
    OPEN(30,FILE="data.dat",STATUS='UNKNOWN',ACTION='READ')
    READ(30,'(t10,i10)') Nx
    READ(30,'(t10,i10)') Nz
    READ(30,'(t10,f10.4)') dx
    READ(30,'(t10,f10.8)') dt
    READ(30,'(t10,f10.4)') fpeak
    READ(30,'(t10,i10)') itmax
    READ(30,'(t10,i10)') lx
    READ(30,'(t10,i10)') lz
    READ(30,'(t10,i10)') nb
    READ(30,'(t10,f10.4)') F
    READ(30,'(t10,i10)') fsf
    READ(30,'(t10,i10)') nphones
    READ(30,'(t10,i10)') npmin_x
    READ(30,'(t10,i10)') npmin_z
    READ(30,'(t10,i10)') dnp_x
    READ(30,'(t10,i10)') dnp_z
    READ(30,'(t10,i10)') nsnaps
    READ(30,'(t10,i10)') snapmin
    READ(30,'(t10,i10)') dsnap
  END SUBROUTINE inputdata


  SUBROUTINE inputmodel(Nx,Nz,C11,C13,C33,C44,rhox)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Nx,Nz
    REAL,    INTENT(OUT), DIMENSION(Nx,Nz) :: C11,C13,C33,C44,rhox
    INTEGER :: i

    OPEN(30,FILE="C11.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx
       READ(30, REC=i) C11(i,:)
    END DO
    CLOSE(30)
    OPEN(30,FILE="C13.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx 
       READ(30, REC=i) C13(i,:)
    END DO
    CLOSE(30)
    OPEN(30,FILE="C33.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx 
       READ(30, REC=i) C33(i,:)
    END DO
    CLOSE(30)
    OPEN(30,FILE="C44.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx 
       READ(30, REC=i) C44(i,:)
    END DO
    CLOSE(30)
    OPEN(30,FILE="rhox.ad",STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='READ', FORM='UNFORMATTED', RECL=Nz*4)
    DO i=1,Nx 
       READ(30, REC=i) rhox(i,:)
    END DO
    CLOSE(30)

  END SUBROUTINE inputmodel


  SUBROUTINE save_shotgather_n_snapshots(csg_ux,csg_uz,csg_P,csg_S,&
       next_ux,next_uz,Sxx,Sxz,Szz,&
       nphones,npmin_x,npmin_z,dnp_x,dnp_z,&
       isnap,nsnaps,snapmin,dsnap,n1,n2,n3,n4,n5,n6,n7,n8,&
       Nx,Nz,it,itmax)
    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(INOUT) :: csg_ux,csg_uz,csg_P,csg_S
    REAL, DIMENSION(:,:), INTENT(IN)    :: next_ux,next_uz,&
         Sxx,Sxz,Szz
    REAL, DIMENSION(Nx) :: aux1
    REAL, DIMENSION(Nz) :: aux2
    INTEGER, INTENT(IN)    :: nphones
    INTEGER, INTENT(IN)    :: npmin_x,npmin_z,dnp_x,dnp_z
    INTEGER, INTENT(INOUT) :: isnap
    INTEGER, INTENT(IN)    :: nsnaps,snapmin,dsnap
    INTEGER, INTENT(IN)    :: n1,n2,n3,n4,n5,n6,n7,n8
    INTEGER, INTENT(IN)    :: Nx,Nz,it,itmax
    INTEGER :: i,ix,iz

    do i=1,nphones
       ix = (i-1)*dnp_x + npmin_x
       iz = (i-1)*dnp_z + npmin_z
       csg_ux(it,i) = next_ux(ix,iz)
       csg_uz(it,i) = next_uz(ix,iz)
       csg_P(it,i) = Sxx(ix,iz) + Szz(ix,iz)
       csg_S(it,i) = Sxz(ix,iz)
    end do

    if (it == itmax) then  
       PRINT*,'Saving CSG'
       do i=1,nphones
          WRITE(n1, REC=i) csg_ux(:,i)
          WRITE(n2, REC=i) csg_uz(:,i)
          WRITE(n5, REC=i) csg_P(:,i)
          WRITE(n7, REC=i) csg_S(:,i)
       end do
    end if

    if (it == (isnap * dsnap) + snapmin) then
       isnap=isnap+1
       if (isnap <= nsnaps) then
          PRINT*, 'Saving snapshot ',isnap,'/',nsnaps
          do i=1,Nx
             WRITE(n3, REC=((isnap -1)*Nx + i)) next_ux(i,:)
             WRITE(n4, REC=((isnap -1)*Nx + i)) next_uz(i,:)
             WRITE(n6, REC=((isnap -1)*Nx + i)) Sxx(i,:)+Szz(i,:)
             WRITE(n8, REC=((isnap -1)*Nx + i)) Sxz(i,:)
          end do
       end if
    end if

  END SUBROUTINE save_shotgather_n_snapshots

END MODULE mdle_io_utils


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



MODULE mdle_source

IMPLICIT NONE

CONTAINS
SUBROUTINE source_init(trs,dt,ns,fpeak,tdelay)
  IMPLICIT NONE
  REAL, INTENT(INOUT), DIMENSION(ns) :: trs
  INTEGER, INTENT(IN) :: ns
  REAL, INTENT(IN) :: dt, fpeak
  REAL, INTENT(OUT) :: tdelay
  INTEGER :: i
  REAL :: t,pi=3.141592653589793238462643383279502884197

  REAL :: wpeak,waux,tt

  wpeak = 2.*pi*fpeak
  waux  = 0.5*wpeak

   tdelay = 6./(5.*fpeak)


  do i=1,ns
     t=(i-1)*dt
     tt = t - tdelay

     trs(i) = exp(-waux*waux*tt*tt/4.)*cos(wpeak*tt)

  end do

END SUBROUTINE source_init

END MODULE
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



