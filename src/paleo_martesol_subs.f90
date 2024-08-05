module paleo_martesol_subs

    use,intrinsic :: iso_fortran_env, only : dp=>real64

    implicit none

    integer :: kmax

    real(dp), parameter :: pi = 4.d0*datan(1.d0)
    real(dp), parameter :: dd = 86400.d0
    real(dp), parameter :: Radm = 3396.19d0
    real(dp), parameter :: uakm = 149597870.7d0
    real(dp), parameter :: epse = 23.43928d0*pi/180.d0
    real(dp), parameter :: a0 = 1.52371034d0
    real(dp), parameter :: GMm = 42828.375816d6
    real(dp), parameter :: GMsol = 1.32712440041279419d20
    real(dp), parameter :: mu = GMm + GMsol
    real(dp), parameter :: norb = dsqrt(mu/(a0**3))
    real(dp), parameter :: Porb = 2.d0*pi/norb
    real(dp), parameter :: Prot = 1.02595676d0*dd
    real(dp), parameter :: fdt = 1.d0/24.d0
    real(dp), parameter :: dt = fdt*Prot


    real(dp) :: t0,intervalo_integracion

    real(dp), allocatable, dimension(:) :: lat,long,h,colat,elong,z,xa,ya,za,wa

    integer, parameter :: mpol = 3
    integer :: nmax
    
    contains

    subroutine calcular_posicion_sol()

        implicit none

        character(len=100) :: fileinDEM,fileinDYN,pathout

        ! Estrategia

        call leer_archivo_entrada(fileinDEM,fileinDYN,pathout)

        call leer_archivo_dem(fileinDEM)

        call leer_archivo_dinamicos(fileinDYN)

        call test()

!        call procesar_datos_escribir_salida(pathout)

    end subroutine

    subroutine test()
        
        character(len=10) :: dia,hora

        integer :: i,j,k,n

        do k=1,1000

            n = 1000 + k

            write(dia,'(I4)') int(n)

            do i = 0,23

                j = 100 + i

                write(hora,'(I3)') int(j)

                print *, 'nombre_de_archivo_'//trim(dia(2:4))//'-'//trim(hora(2:3))//'.out'
            
            end do

        end do

        
    end subroutine test

    subroutine procesar_datos_escribir_salida(pathout)

        implicit none

        real(dp) :: t,e,eps,pibar
        real(dp) :: xeq,yeq,zeq
        real(dp) :: Ac,Alt,dSol
        integer :: k,d,hh

        character(len=100), intent(in) :: pathout
        character(len=100) ::  fileout

        call buscar_e_y_eps(t0,e,eps,pibar)

        t = 0
        d = 1
        hh = 0

        do while (t.le.Porb)

            fileout = trim(pathout)//'.out'
        
            open(unit=12,file=trim(pathout))

            call calcular_coords_eqs_sol(t,e,eps,pibar,xeq,yeq,zeq)

            do k=1,size(colat)

                call calcular_Acimut_Alt_sol(k,xeq,yeq,zeq,Ac,Alt,dSol)

                write(12,*) long(k),lat(k),Ac,Alt,dSol

            end do

            t = t + dt
            h = h + 1



        end do
        
    end subroutine procesar_datos_escribir_salida

    subroutine calcular_coords_eqs_sol(t,e,eps,pibar,xeq,yeq,zeq)

        real(dp), intent(in) :: t,e,eps,pibar

        real(dp), intent(out) :: xeq,yeq,zeq

        real(dp) :: AM,AE,rr,xo,yo
        real(dp) :: xecl,yecl,zecl

        AM = norb*t

        call SOLKEP(e,AM,AE)

        rr =  a0*(1.D0-e*dcos(AE))

        xo = a0*(dcos(AE) - e)
        yo = a0*dsqrt(1.d0-e*e)*dsin(AE)

        xecl = xo*dcos(pibar) - yo*dsin(pibar)
        yecl = xo*dsin(pibar) + yo*dcos(pibar)
        zecl = 0.d0

        xeq = xecl
        yeq = yecl*dcos(eps) - zecl*dsin(eps)
        zeq = yecl*dsin(eps) + zecl*dcos(eps)

    end subroutine

    subroutine calcular_Acimut_Alt_sol(k,xeq,yeq,zeq,Ac,Alt,dSol)

        integer, intent(in) :: k

        real(dp), intent(in) :: xeq,yeq,zeq
        real(dp), intent(out) :: Ac,Alt,dSol

        real(dp) :: xsup,ysup,zsup
        real(dp) :: xsol,ysol,zsol

    
        xsup = (Radm + z(k))*dsin(colat(k))*dcos(elong(k))
        ysup = (Radm + z(k))*dsin(colat(k))*dsin(elong(k))
        zsup = (Radm + z(k))*dcos(colat(k))

        xsol = - (xsup + xeq)
        ysol = - (ysup + yeq)
        zsol = - (zsup + zeq)

        dSol = PYTHAG(PYTHAG(xsol,ysol),zsol)

        Ac = datan2(xsol,ysol)
        Ac = Ac*180.d0/pi
        if (Ac.lt.0.d0) Ac = Ac + 360.d0
        Alt = datan2(zsol,pythag(xsol,ysol))
        Alt = Alt*180.d0/pi
        
    end subroutine calcular_Acimut_Alt_sol

    subroutine leer_archivo_entrada(fileinDEM,fileinDYN,pathout)

        implicit none

        character(len=100) :: algo
        character(len=100),intent(out) :: fileinDEM,fileinDYN,pathout
        
        open(unit=10,file="../in/PALEOMARTESOL.IN",status="unknown")

        read(10,*) algo
        read(10,*) algo
        read(10,*) t0
        read(10,*) algo
        read(10,*) fileinDYN
        read(10,*) algo
        read(10,*) intervalo_integracion
        read(10,*) algo
        read(10,*) fileinDEM
        read(10,*) algo
        read(10,*) pathout
        
    end subroutine leer_archivo_entrada

    subroutine leer_archivo_dem(filein)
        
        implicit none

        integer :: k,feof

        real(dp) :: xdata,ydata,zdata

        character(len=100),intent(in) :: filein

        open(unit=11,file=trim(filein))

        k = 1
        
        feof = 0
        
        print *, 'Procesando archivo "'//trim(filein)//'"'

        do while(feof.eq.0)

            read(11,*,iostat=feof) xdata,ydata,zdata

            if (feof.gt.0) then
                print *,'Revisar archivo de entrada'
                exit
            else if (feof.lt.0) then
                print *,  '¡Listo!'
                exit
            end if
        
            k = k + 1
        
        end do 

        close(11)
        
        kmax = k - 1

        allocate(colat(kmax), source=0.d0)
        allocate(elong(kmax), source=0.d0)
        allocate(z(kmax), source=0.d0)

        allocate(long(kmax), source=0.d0)
        allocate(lat(kmax), source=0.d0)
        allocate(h(kmax), source=0.d0)

        open(unit=11,file=trim(filein))

        do k=1,kmax

            read(11,*,iostat=feof) long(k),lat(k),h(k)

            colat(k) = 90.d0 - lat(k)

            if (long(k).ge.0.d0) then
                elong(k) = long(k)
            else
                elong(k) = 360.d0 + long(k)
            end if

            colat(k) = colat(k)*pi/180.d0
            elong(k) = elong(k)*pi/180.d0
            z(k) = h(k)*1.d-3

        end do
        
    end subroutine leer_archivo_dem

    subroutine leer_archivo_dinamicos(filein)
        
        implicit none

        integer :: k,feof

        real(dp) :: t,e,eps,pibar

        character(len=100),intent(in) :: filein

        open(unit=11,file=trim(filein))

        k = 1
        
        feof = 0
        
        print *, 'Procesando archivo "'//trim(filein)//'"'

        do while(feof.eq.0)

            read(11,*,iostat=feof) t,e,eps,pibar

            if (feof.gt.0) then
                print *,'Revisar archivo de entrada'
                exit
            else if (feof.lt.0) then
                print *,  '¡Listo!'
                exit
            end if
        
            k = k + 1
        
        end do 

        close(11)
        
        kmax = k - 1

        nmax = kmax

        allocate(xa(kmax), source=0.d0)
        allocate(ya(kmax), source=0.d0)
        allocate(za(kmax), source=0.d0)
        allocate(wa(kmax), source=0.d0)

        open(unit=11,file=trim(filein))

        do k=1,kmax

            read(11,*,iostat=feof) t,e,eps,pibar

            xa(k) = t
            ya(k) = e
            za(k) = eps
            wa(k) = pibar

        end do
        
    end subroutine leer_archivo_dinamicos

    subroutine buscar_e_y_eps(t_ini,e,eps,pibar)

        real(dp), intent(in) :: t_ini
        real(dp),intent(out) :: e,eps,pibar
        real(dp) :: y,dy
        integer :: j,k

        call hunt(xa,nmax,t_ini,j)
        k = min(max(j-(mpol-1)/2,1),nmax+1-mpol)

        call polint(xa(k),ya(k),mpol,t_ini,y,dy)

        e = y

        call polint(xa(k),za(k),mpol,t_ini,y,dy)

        eps = y

        call polint(xa(k),wa(k),mpol,t_ini,y,dy)

        pibar = y
        
    end subroutine buscar_e_y_eps

!=======================================================================
    SUBROUTINE hunt(xx,n,x,jlo)
        INTEGER :: jlo,n
        REAL (kind=dp) x,xx(n)
        INTEGER :: inc,jhi,jm
        LOGICAL :: ascnd
        ascnd=xx(n).gt.xx(1)
        if(jlo.le.0.or.jlo.gt.n)then
          jlo=0
          jhi=n+1
          goto 3
        endif
        inc=1
        if(x.ge.xx(jlo).eqv.ascnd)then
1         jhi=jlo+inc
          if(jhi.gt.n)then
            jhi=n+1
          else if(x.ge.xx(jhi).eqv.ascnd)then
            jlo=jhi
            inc=inc+inc
            goto 1
          endif
        else
          jhi=jlo
2         jlo=jhi-inc
          if(jlo.lt.1)then
            jlo=0
          else if(x.lt.xx(jlo).eqv.ascnd)then
            jhi=jlo
            inc=inc+inc
            goto 2
          endif
        endif
3       if(jhi-jlo.eq.1)return
        jm=(jhi+jlo)/2
        if(x.gt.xx(jm).eqv.ascnd)then
          jlo=jm
        else
          jhi=jm
        endif
        goto 3
    END SUBROUTINE hunt

    subroutine polint(xp,yp,n,x,y,dy)
        integer :: n,k,m,ns
        real (kind=dp), intent(in) :: x,xp(n),yp(n)
        real (kind=dp), intent(out) :: y,dy
        integer, parameter :: kmaxi=10
        real (kind=dp) :: dlen,dif,dift,ho,hp,w,c(kmaxi),d(kmaxi)
        ns=1
        dif=dabs(x-xp(1))
        do k=1,n
            dift=dabs(x-xp(k))
            if (dift.lt.dif) then
            ns=k
            dif=dift
            endif
            c(k)=yp(k)
            d(k)=yp(k)
        end do
        y=yp(ns)
        ns=ns-1
        do m=1,n-1
            do k=1,n-m
            ho=xp(k)-x
            hp=xp(k+m)-x
            w=c(k+1)-d(k)
            dlen=ho-hp
            if(dlen.eq.0.d0) print *,'failure in polint'
            dlen=w/dlen
            d(k)=hp*dlen
            c(k)=ho*dlen
            end do
            if (2*ns.lt.n-m)then
            dy=c(ns+1)
            else
            dy=d(ns)
            ns=ns-1
            endif
            y=y+dy
        end do
        return
    end subroutine polint

    FUNCTION PYTHAG(A,B)
        REAL*8 A,B,PYTHAG
        REAL*8 ABSA,ABSB
        ABSA=ABS(A)
        ABSB=ABS(B)
        IF(ABSA.GT.ABSB)THEN
            PYTHAG=ABSA*SQRT(1.D0+(ABSB/ABSA)**2)
        ELSE
            IF(ABSB.EQ.0.D0)THEN
            PYTHAG=0.D0
            ELSE
            PYTHAG=ABSB*SQRT(1.D0+(ABSA/ABSB)**2)
            ENDIF
        ENDIF
        RETURN
    END

! SOLUCION DE LA ECUACION DE KEPLER ELIPTICA:
    SUBROUTINE SOLKEP(EX,M,E)
!
! SOLUCION ITERATIVA DE LA ECUACION DE KEPLER
! ENTRA:EX   EXCENTRICIDAD            (<1)
!       M    ANOMALIA MEDIA           (RADIANES)
! SALE: E    ANOMALIA EXCENTRICA      (RADIANES)
!
    IMPLICIT REAL*8 (A-H,O-Z)
    REAL*8 M,MK
    integer :: NITER,NDIC

    TOLE=1.D-10
    DPI=8.D0*DATAN(1.D0)
    M=DMOD(M,DPI)
    E=M
    NITER=0

100     E0=E
    
    SE=DSIN(E0)
    CE=DCOS(E0)
    ES=EX*SE
    EC=1.D0-EX*CE
    MK=E0-ES
    U=(MK-M)/EC
    XPRI=E0-U
    XSEG=E0-U/(1.D0-U*ES)
    E=(XPRI+XSEG)/2.D0
    DEX=DABS(E-E0)
    
    NITER=NITER+1
    IF(NITER.GT.20)GOTO 200
    IF(DEX.GT.TOLE)GOTO 100                
    RETURN

! SI EL NUMERO DE ITERACIONES ES > 20 PRUEBA CON BISECCION:
    
!   METODO DICOTOMICO:
200      CONTINUE        
    NDIC=0
    E0=-DPI
    DE0=DPI/10.D0

400     DE=DE0/(10.D0**NDIC)
    SE=DSIN(E0)
    CE=DCOS(E0)
    ES=EX*SE
    EM0=E0-ES-M
    NITER=0
    
300      E1=E0+DE
    NITER=NITER+1
    
    IF(NITER.GT.100)THEN
!      WRITE(50,*)'ERROR EN LA SOLUCION DE LA ECUACION DE KEPLER'
    RETURN
    ENDIF

    SE=DSIN(E1)
    CE=DCOS(E1)
    ES=EX*SE
    EM1=E1-ES-M
    IF(EM1*EM0.GT.0.D0)THEN
    E0=E1
    EM0=EM1
    GOTO 300
    ELSE
    NDIC=NDIC+1
    IF(NDIC.EQ.3)THEN
    E=E1
    RETURN
    ENDIF
    GOTO 400
    
    ENDIF

    RETURN
    END
    
end module paleo_martesol_subs