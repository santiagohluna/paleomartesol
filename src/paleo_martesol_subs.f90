module paleo_martesol_subs

    use,intrinsic :: iso_fortran_env, only : dp=>real64

    implicit none

    integer :: kmax

    real(dp), parameter :: pi = 4.d0*datan(1.d0)
    real(dp), parameter :: dd = 86400.d0
    real(dp), parameter :: Radm = 3396.19d0
    real(dp), parameter :: uakm = 149597870.7d0
    real(dp), parameter :: epse = 23.43928d0*pi/180.d0
    real(dp), parameter :: a0ua = 1.52371034d0
    real(dp), parameter :: a0km = a0ua*uakm
    real(dp), parameter :: a0 = a0km*1.d3
    real(dp), parameter :: GMm = 42828.375816d9
    real(dp), parameter :: GMsol = 1.32712440041279419d20
    real(dp), parameter :: mu = GMm + GMsol
    real(dp), parameter :: norb = dsqrt(mu/(a0**3))
    real(dp), parameter :: Porb = 2.d0*pi/norb
    real(dp), parameter :: Prot = 1.02595676d0*dd
    real(dp), parameter :: thp = 2.d0*pi/Prot
    real(dp), parameter :: fdt = 1.d0/24.d0
    real(dp), parameter :: dt = fdt*Prot

    character, parameter :: cr = achar(13)

    real(dp) :: t0,intervalo_integracion,tfin,theta0

    real(dp), allocatable, dimension(:) :: lat,long,h,colat,elong,z,xa,ya,za,wa

    integer, parameter :: mpol = 3
    integer :: nmax
    
    contains

    subroutine calcular_posicion_sol()

        implicit none

        integer :: modelo_dinamico
        character(len=100) :: fileinDEM,pathout

        ! Estrategia

        call leer_archivo_entrada(fileinDEM,modelo_dinamico,pathout)

        call leer_archivo_dem(fileinDEM)

        call leer_archivo_dinamicos(modelo_dinamico)

        call procesar_o_probar(pathout)

    end subroutine

    subroutine procesar_o_probar(pathout)

        implicit none

        character(len=1) :: entrada

        character(len=100) :: pathout

        logical :: ltest,lproc

        print *,'Ingrese P para procesar o T para probar:'
        
        read *, entrada

        lproc = (entrada.eq.'P').or.(entrada.eq.'p')
        ltest = (entrada.eq.'T').or.(entrada.eq.'t')

        if(lproc) then

            call procesar_datos_escribir_salida(pathout)

        else if(ltest) then

            call test()

        else

            print *,'Debe ingresar la letra P/p o T/t'

        end if

        
    end subroutine procesar_o_probar

    subroutine test()
        
        character(len=10) :: dia,hora,anio

        integer :: d,hh,yy,counter,max
        real(dp) :: t

        t = 0
        d = 1
        hh = 0
        yy = 1

        counter = 0
        max = int(intervalo_integracion*24*Porb/Prot)

        print *,'max = ',max

        tfin = intervalo_integracion*Porb

        do while (t <= tfin)

            write(dia,'(I4)') 1000 + d
            
            write(hora,'(I3)') 100 + hh

            write(anio,'(I3)') 100 + yy

            if (intervalo_integracion > 1) then

                print *,trim(anio(2:3))//'-'//trim(dia(2:4))//'-'//trim(hora(2:3))//'.out'
            
            else
    
                print *,trim(dia(2:4))//'-'//trim(hora(2:3))//'.out'
    
            end if

            t = t + dt
            hh = hh + 1

            if(hh > 23) then 

                d = d + 1 
                hh = 0

            end if

            if (t > yy*Porb) then
                 
                yy = yy + 1
                d = 1
                hh = 0

            end if

            call progress_bar(counter,max)

            counter = counter + 1

        end do

        print *,'counter = ',counter

        
    end subroutine test

    subroutine procesar_datos_escribir_salida(pathout)

        implicit none

        integer :: counter,max

        real(dp) :: t,e,eps,pibar
        real(dp) :: xeq,yeq,zeq,Ls

        integer :: d,hh,yy

        character(len=100), intent(in) :: pathout

        counter = 0
        max = int(intervalo_integracion*24*Porb/Prot)        

        call buscar_e_y_eps(t0,e,eps,pibar)

        t = 0
        d = 1
        hh = 0
        yy = 1

        tfin = intervalo_integracion*Porb

        do while (t.le.tfin)

            call calcular_coords_eqs_sol(t,e,eps,pibar,xeq,yeq,zeq,Ls)

            call escribir_archivo_salida(pathout,xeq,yeq,zeq,Ls,d,hh,yy)

            call progress_bar(counter,max)

            counter = counter + 1

            t = t + dt
            hh = hh + 1

            if(hh > 23) then 

                d = d + 1 
                hh = 0

            end if

            if (t > yy*Porb) then
                 
                yy = yy + 1
                d = 1
                hh = 0

            end if

        end do

        print *,' '
        print *,' '
        print *,'¡Listo!'
        print *,' '
        
    end subroutine procesar_datos_escribir_salida

    subroutine escribir_archivo_salida(pathout,xeq,yeq,zeq,Ls,d,hh,yy)

        implicit none

        character(len=100), intent(in) :: pathout
        real(dp), intent(in) :: xeq,yeq,zeq,Ls
        integer, intent(in) :: d,hh,yy

        integer :: k
        character(len=10) :: dia,hora,anio
        character(len=100) :: fileout
        real(dp) :: Ac,Alt,dSol

        write(dia,'(I4)') 1000 + d
            
        write(hora,'(I3)') 100 + hh
        
        write(anio,'(I3)') 100 + yy

        if (intervalo_integracion > 1) then

            fileout = trim(pathout)//trim(anio(2:3))//'-'//trim(dia(2:4))//'-'//trim(hora(2:3))//'.out'
        
        else

            fileout = trim(pathout)//trim(dia(2:4))//'-'//trim(hora(2:3))//'.out'

        end if
    
        open(unit=12,file=fileout)
        open(unit=13,file='../out/analema.out',status='unknown')

        do k=1,kmax

            call calcular_Acimut_Alt_sol(k,xeq,yeq,zeq,Ac,Alt,dSol)

            write(12,*) long(k),lat(k),Ac,Alt,dSol

            if((k==1).and.(hh==12)) then 

                if (intervalo_integracion > 1) then

                    write(13,*) Ac,Alt,Ls,d,yy
                
                else
        
                    write(13,*) Ac,Alt,Ls,d
        
                end if

            end if

        end do
        
    end subroutine escribir_archivo_salida

    subroutine calcular_coords_eqs_sol(t,e,eps,pibar,xeq,yeq,zeq,Ls)

        real(dp), intent(in) :: t,e,eps,pibar

        real(dp), intent(out) :: xeq,yeq,zeq,Ls

        real(dp) :: AM,AE,xo,yo,theta
        real(dp) :: xecl,yecl,zecl
        real(dp) :: xeq_cel,yeq_cel,zeq_cel

        AM = PRIVUE(norb*t)
        theta = theta0 + thp*t

        call SOLKEP(e,AM,AE)

        xo = a0ua*(dcos(AE) - e)
        yo = a0ua*dsqrt(1.d0-e*e)*dsin(AE)

        xecl = xo*dcos(pibar) - yo*dsin(pibar)
        yecl = xo*dsin(pibar) + yo*dcos(pibar)
        zecl = 0.d0

        Ls = atg(-yecl,-xecl)*180.d0/pi ! Cálculo de la longitud solar

        xeq_cel = xecl
        yeq_cel = yecl*dcos(eps) - zecl*dsin(eps)
        zeq_cel = yecl*dsin(eps) + zecl*dcos(eps)

        xeq = xeq_cel*dcos(theta) + yeq_cel*dsin(theta)
        yeq = -xeq_cel*dsin(theta) + yeq_cel*dcos(theta)
        zeq = zeq_cel

    end subroutine

    subroutine calcular_Acimut_Alt_sol(k,xeq,yeq,zeq,Ac,Alt,dSol)

        integer, intent(in) :: k

        real(dp), intent(in) :: xeq,yeq,zeq
        real(dp), intent(out) :: Ac,Alt,dSol

        real(dp) :: xsup,ysup,zsup
        real(dp) :: xsol_eq,ysol_eq,zsol_eq
        real(dp) :: xsol_h,ysol_h,zsol_h,xsol_h1,ysol_h1,zsol_h1,xsol_h2,ysol_h2,zsol_h2

    
        xsup = (Radm + z(k))*dsin(colat(k))*dcos(elong(k))/uakm
        ysup = (Radm + z(k))*dsin(colat(k))*dsin(elong(k))/uakm
        zsup = (Radm + z(k))*dcos(colat(k))/uakm

        xsol_eq = - (xsup + xeq)
        ysol_eq = - (ysup + yeq)
        zsol_eq = - (zsup + zeq)

        dSol = PYTHAG(PYTHAG(xsol_eq,ysol_eq),zsol_eq)

        xsol_h1 =  xsol_eq*dcos(elong(k)) + ysol_eq*dsin(elong(k))
        ysol_h1 = -xsol_eq*dsin(elong(k)) + ysol_eq*dcos(elong(k))
        zsol_h1 =  zsol_eq

        xsol_h2 = xsol_h1*dcos(colat(k)) - zsol_h1*dsin(colat(k))
        ysol_h2 = ysol_h1
        zsol_h2 = xsol_h1*dsin(colat(k)) + zsol_h1*dcos(colat(k))

        xsol_h = -xsol_h2
        ysol_h = -ysol_h2
        zsol_h =  zsol_h2

        Ac = 360.d0 - atg(ysol_h,xsol_h)*180.d0/pi
        Alt = datan2(zsol_h,pythag(xsol_h,ysol_h))
        Alt = Alt*180.d0/pi
        
    end subroutine calcular_Acimut_Alt_sol

    subroutine leer_archivo_entrada(fileinDEM,modelo_dinamico,pathout)

        implicit none

        integer, intent(out) :: modelo_dinamico
        character(len=100) :: algo
        character(len=100),intent(out) :: fileinDEM,pathout
        
        open(unit=10,file="../in/PALEOMARTESOL.IN",status="unknown")

        read(10,*) algo
        read(10,*) algo
        read(10,*) t0
        read(10,*) algo
        read(10,*) theta0
        read(10,*) algo
        read(10,*) algo
        read(10,*) modelo_dinamico
        read(10,*) algo
        read(10,*) intervalo_integracion
        read(10,*) algo
        read(10,*) fileinDEM
        read(10,*) algo
        read(10,*) pathout

        if(intervalo_integracion>19) then 
            print *,'El intervalo de integración debe ser menor o igual a 99 periodos orbitales de Marte.'
            call exit()
        end if
        
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

    subroutine leer_archivo_dinamicos(id_dinamico)
        
        implicit none

        integer, intent(in) :: id_dinamico

        integer :: n,feof

        real(dp) :: t,e,eps,pibar

        character(len=100):: filein

        if(t0 < 0) then 
            if (id_dinamico == 3) then
                filein = '../data/INSOLN.LA2003.MARS.ASC'
            else if (id_dinamico == 4) then
                filein = '../data/INSOLN.LA2004.MARS.ASC'
            else
                print *,'Error en el identificador del modelo dinámico.'
            end if
        else
            if (id_dinamico == 3) then
                filein = '../data/INSOLP.LA2003.MARS.ASC'
            else if (id_dinamico == 4) then 
                filein = '../data/INSOLP.LA2004.MARS.ASC'
            else
                print *,'Error en el identificador del modelo dinámico.'
            end if
        end if

        open(unit=11,file=trim(filein))

        n = 1
        
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
        
            n = n + 1
        
        end do 

        close(11)
        
        nmax = n - 1 

        allocate(xa(nmax), source=0.d0)
        allocate(ya(nmax), source=0.d0)
        allocate(za(nmax), source=0.d0)
        allocate(wa(nmax), source=0.d0)

        open(unit=11,file=trim(filein))

        do n=1,nmax

            read(11,*,iostat=feof) t,e,eps,pibar

            xa(n) = t
            ya(n) = e
            za(n) = eps
            wa(n) = pibar

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

    real(dp) function atg(y,x) result(retval)

        real(dp), intent(in) :: x,y

        retval = atan(y,x)

        if (y < 0.0d0) retval = retval + 2.d0*pi
        
    end function atg

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

    FUNCTION PRIVUE(X)

      IMPLICIT NONE

      REAL*8 X,PRIVUE

      IF (X < 0.D0) X = X + 2.d0*pi
      IF (X > 2.d0*pi) THEN
       PRIVUE = X - DBLE(INT(X/2.d0*pi))*2.d0*pi
      ELSE
       PRIVUE = X
      END IF

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
    DPI=2.d0*pi
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

    subroutine progress_bar(iteration, maximum)
    !
    ! Prints progress bar.
    !
    ! Args: 
    !     iteration - iteration number
    !     maximum - total iterations
    !
        implicit none
        integer :: iteration, maximum
        integer :: counter
        integer :: step, done

        step = nint(iteration * 100 / (1.0 * maximum))
        done = floor(step / 10.0)  ! mark every 10%

        do counter = 1, 36                    ! clear whole line - 36 chars
            write(6,'(a)',advance='no') '\b'  ! (\b - backslash)
        end do

        write(6,'(a)',advance='no') ' -> In progress... ['
        if (done .LE. 0) then
            do counter = 1, 10
                write(6,'(a)',advance='no') '='
            end do
        else if ((done .GT. 0) .and. (done .LT. 10)) then
            do counter = 1, done
                write(6,'(a)',advance='no') '#'
            end do
            do counter = done+1, 10
                write(6,'(a)',advance='no') '='
            end do 
        else
            do counter = 1, 10
                write(6,'(a)',advance='no') '#'
            end do
        end if
        write(6,'(a)',advance='no') '] '
        write(6,'(I3.1)',advance='no') step
        write(6,'(a)',advance='no') '%'
    end
    
end module paleo_martesol_subs