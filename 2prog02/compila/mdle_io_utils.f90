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


