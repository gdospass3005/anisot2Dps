MODULE mdle_source

  IMPLICIT NONE

CONTAINS

  SUBROUTINE source_dervgauss_init(trs,dt,ns,fpeak,tdelay)
    IMPLICIT NONE
    REAL, INTENT(INOUT), DIMENSION(ns) :: trs
    INTEGER, INTENT(IN) :: ns
    REAL, INTENT(IN) :: dt, fpeak
    REAL, INTENT(OUT) :: tdelay
    INTEGER :: i
    REAL :: bw,t,pi=3.1415926536

    tdelay = 4.0 /(3.*fpeak*SQRT(2.*pi))

    bw=1.
    do i=1,ns
       t=(i-1)*dt - tdelay
       trs(i) = -dervgauss(t,fpeak,bw)
    end do

    OPEN(81,FILE='sourcememory.ad',STATUS='UNKNOWN',ACCESS='DIRECT', &
         ACTION='WRITE', FORM='UNFORMATTED', RECL=Ns*4)
    WRITE(81, REC=1) trs(:)
    CLOSE(81)

  END SUBROUTINE source_dervgauss_init


  REAL FUNCTION dervgauss(t,fpeak,bombweight)
    REAL, INTENT(IN) :: t,fpeak
    REAL :: x,xx,pi=3.1415926536
    REAL, INTENT(IN) :: bombweight

    x=3*fpeak*t
    xx=x*x

    dervgauss= bombweight*(-2.*pi*t)*EXP(-pi*xx)

  END FUNCTION dervgauss


  SUBROUTINE source_UU3_2Dp(t,it,trs,s,r,tdelay)
    IMPLICIT NONE
    REAL, INTENT(IN) :: t, tdelay
    INTEGER, INTENT(IN) :: it
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT), DIMENSION(:,:) :: s
    REAL, INTENT(IN), DIMENSION(:) :: trs
    REAL :: stemp,lx,lz
    INTEGER :: ix,iz,ir,ixlow,ixhigh,izlow,izhigh
    REAL :: ttot,tmin
    REAL :: PI=3.1415926536,r2,d2,wd,aux,dist

    s=0.

    IF (t <= (2.0*tdelay)) THEN
       ir = CEILING(r)
       lx = ir - 0.5
       lz = ir - 0.5
       stemp = trs(it)
       IF (r /= 0) THEN
          ixlow =  1
          ixhigh = 2*ir + 1
          izlow =  1
          izhigh = 2*ir + 1
          r2 = (REAL(r))**2
          aux = LOG(0.25) / r2 
          DO ix=ixlow,ixhigh
             DO iz=izlow,izhigh
                d2 = (( REAL(lx-ix) )**2 + ( REAL(lz-iz) )**2)
                dist=sqrt(d2)
                IF (d2 <= r2) THEN
                   wd = EXP( d2 * aux  )
                   !wd = (r - dist)/r
!                   IF (dist /= 0) THEN
                      ! s(ix,iz) = ((iz-lz)/dist)*stemp * wd*dist
                      ! s(ix,iz) = ((iz-lz)/dist)*stemp * wd
                      s(ix,iz) = stemp * wd
!                   ELSE
!                      s(ix,iz) = 0
!                   END IF
                END IF
             END DO
          END DO
       ELSE
          STOP 'source radius must be greater or equal to 1.'
       END IF
    END IF

  END SUBROUTINE source_UU3_2Dp



  SUBROUTINE source_UU1_2Dc(t,it,trs,s,r,tdelay)
    IMPLICIT NONE
    REAL, INTENT(IN) :: t, tdelay
    INTEGER, INTENT(IN) :: it
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT), DIMENSION(:,:) :: s
    REAL, INTENT(IN), DIMENSION(:) :: trs
    REAL :: stemp,lx,lz
    INTEGER :: ix,iz,ir,ixlow,ixhigh,izlow,izhigh
    REAL :: ttot,tmin
    REAL :: PI=3.1415926536,r2,d2,wd,aux,dist

    s=0.

    IF (t <= (2.0*tdelay)) THEN
       ir = CEILING(r)
       lx = ir + 0.5
       lz = ir + 0.5
       stemp = trs(it)
       IF (r /= 0) THEN
          ixlow =  1
          ixhigh = 2*ir + 1
          izlow =  1
          izhigh = 2*ir + 1
          r2 = (REAL(r))**2
          aux = LOG(0.25) / r2 
          DO ix=ixlow,ixhigh
             DO iz=izlow,izhigh
                d2 = (( REAL(lx-ix) )**2 + ( REAL(lz-iz) )**2)
                dist=sqrt(d2)
                IF (d2 < r2) THEN
                   wd = EXP( d2 * aux  )
                   ! wd = (r - dist)/r
                   IF (dist /= 0) THEN
                      s(ix,iz) = ((ix-lx)/dist)*stemp*wd*dist
                   ELSE
                      s(ix,iz) = 0
                   END IF
                END IF
             END DO
          END DO
       ELSE
          STOP 'source radius must be greater or equal to 1.'
       END IF
    END IF


  END SUBROUTINE source_UU1_2Dc



  SUBROUTINE source_UU3_2Dc(t,it,trs,s,r,tdelay)
    IMPLICIT NONE
    REAL, INTENT(IN) :: t, tdelay
    INTEGER, INTENT(IN) :: it
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT), DIMENSION(:,:) :: s
    REAL, INTENT(IN), DIMENSION(:) :: trs
    REAL :: stemp,lx,lz
    INTEGER :: ix,iz,ir,ixlow,ixhigh,izlow,izhigh
    REAL :: ttot,tmin
    REAL :: PI=3.1415926536,r2,d2,wd,aux,dist

    s=0.

    IF (t <= (2.0*tdelay)) THEN
       ir = CEILING(r)
       ! lx = ir - 0.5
       ! lz = ir - 0.5
       lx = ir + 0.5
       lz = ir + 0.5
       stemp = trs(it)
       IF (r /= 0) THEN
          ixlow =  1
          ixhigh = 2*ir + 1
          izlow =  1
          izhigh = 2*ir + 1
          r2 = (REAL(r))**2
          aux = LOG(0.25) / r2 
          DO ix=ixlow,ixhigh
             DO iz=izlow,izhigh
                d2 = (( REAL(lx-ix) )**2 + ( REAL(lz-iz) )**2)
                dist=sqrt(d2)
                IF (d2 <= r2) THEN
                   wd = EXP( d2 * aux  )
                   !wd = (r - dist)/r
                   IF (dist /= 0) THEN
                      s(ix,iz) = ((iz-lz)/dist)*stemp * wd*dist
                   ELSE
                      s(ix,iz) = 0
                   END IF
                END IF
             END DO
          END DO
       ELSE
          STOP 'source radius must be greater or equal to 1.'
       END IF
    END IF

  END SUBROUTINE source_UU3_2Dc



  SUBROUTINE source_UU1_2D(t,it,trs,s,r,tdelay)
    IMPLICIT NONE
    REAL, INTENT(IN) :: t, tdelay
    INTEGER, INTENT(IN) :: it
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT), DIMENSION(:,:) :: s
    REAL, INTENT(IN), DIMENSION(:) :: trs
    REAL :: stemp,lx,lz
    INTEGER :: ix,iz,ir,ixlow,ixhigh,izlow,izhigh
    REAL :: ttot,tmin
    REAL :: PI=3.1415926536,r2,d2,wd,aux,dist

    s=0.

    IF (t <= (2.0*tdelay)) THEN
       ir = CEILING(r)
       lx = ir + 0.
       lz = ir + 0.
       stemp = trs(it)

       IF (r /= 0) THEN
          ixlow =  1
          ixhigh = 2*ir + 1
          izlow =  1
          izhigh = 2*ir + 1
          r2 = (REAL(r))**2
          aux = LOG(0.25) / r2 
          DO ix=ixlow,ixhigh
             DO iz=izlow,izhigh
                d2 = (( REAL(lx-ix) )**2 + ( REAL(lz-iz) )**2)
                dist=sqrt(d2)
                IF (d2 < r2) THEN
                   wd = EXP( d2 * aux  )
                   ! wd = (r - dist)/r
                   IF (dist /= 0) THEN
                      s(ix,iz) = ((ix-lx)/dist)*stemp*wd*dist
                   ELSE
                      s(ix,iz) = 0
                   END IF
                END IF
             END DO
          END DO
       ELSE
          STOP 'source radius must be greater or equal to 1.'
       END IF
    END IF

  END SUBROUTINE source_UU1_2D



  SUBROUTINE source_UU3_2D(t,it,trs,s,r,tdelay)
    IMPLICIT NONE
    REAL, INTENT(IN) :: t, tdelay
    INTEGER, INTENT(IN) :: it
    REAL, INTENT(IN) :: r
    REAL, INTENT(OUT), DIMENSION(:,:) :: s
    REAL, INTENT(IN), DIMENSION(:) :: trs
    REAL :: stemp,lx,lz
    INTEGER :: ix,iz,ir,ixlow,ixhigh,izlow,izhigh
    REAL :: ttot,tmin
    REAL :: PI=3.1415926536,r2,d2,wd,aux,dist

    s=0.

    IF (t <= (2.0*tdelay)) THEN
       ir = CEILING(r)
       lx = ir - 0.5
       lz = ir - 0.5
       stemp = trs(it)
       IF (r /= 0) THEN
          ixlow =  1
          ixhigh = 2*ir + 1
          izlow =  1
          izhigh = 2*ir + 1
          r2 = (REAL(r))**2
          aux = LOG(0.25) / r2 
          DO ix=ixlow,ixhigh
             DO iz=izlow,izhigh
                d2 = (( REAL(lx-ix) )**2 + ( REAL(lz-iz) )**2)
                dist=sqrt(d2)
                IF (d2 <= r2) THEN
                   wd = EXP( d2 * aux  )
                   !wd = (r - dist)/r
                   IF (dist /= 0) THEN
                      s(ix,iz) = ((iz-lz)/dist)*stemp * wd*dist
                   ELSE
                      s(ix,iz) = 0
                   END IF
                END IF
             END DO
          END DO
       ELSE
          STOP 'source radius must be greater or equal to 1.'
       END IF
    END IF

  END SUBROUTINE source_UU3_2D





END MODULE mdle_source

