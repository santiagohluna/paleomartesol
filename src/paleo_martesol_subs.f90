module paleo_martesol_subs

    use,intrinsic :: iso_fortran_env, only : dp=>real64

    implicit none

    integer :: kmax

    real(dp), parameter :: pi = 4.d0*datan(1.d0)
    real(dp), parameter :: radm = 3396.19d0
    real(dp), parameter :: uakm = 149597870.7d0
    real(dp), parameter :: epse = 23.43928d0*pi/180.d0

    real(dp) :: t0,dt,fdt

    real(dp), allocatable, dimension(:) :: colat,elong,z
    
    contains

    subroutine calcular_posicion_sol()

        implicit none

        character(len=100) :: filein,pathout

        ! Estrategia

        call leer_archivo_entrada(filein,pathout)

        call leer_archivo_dem(filein)

        ! call procesar_datos_escribir_salida(pathout)

    end subroutine

    subroutine leer_archivo_entrada(filein,pathout)

        implicit none

        character(len=100) :: algo
        character(len=100),intent(out) :: filein,pathout
        
        open(unit=10,file="../in/PALEOMARTESOL.IN",status="unknown")

        read(10,*) algo
        read(10,*) algo
        read(10,*) t0
        read(10,*) algo
        read(10,*) dt
        read(10,*) algo
        read(10,*) fdt
        read(10,*) algo
        read(10,*) filein
        read(10,*) algo
        read(10,*) pathout
        
    end subroutine leer_archivo_entrada

    subroutine leer_archivo_dem(filein)
        
        implicit none

        integer :: k,feof

        real(dp) :: long,lat,h

        character(len=100),intent(in) :: filein

        open(unit=11,file=trim(filein))

        k = 1
        
        feof = 0
        
        print *, 'Procesando archivo "'//trim(filein)//'"'

        do while(feof.eq.0)

            read(11,*,iostat=feof) long,lat,h

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

        open(unit=11,file=trim(filein))

        do k=1,kmax

            read(11,*,iostat=feof) long,lat,h

            colat(k) = 90.d0 - lat

            if (long.ge.0.d0) then
                elong(k) = long
            else
                elong(k) = 360.d0 + long
            end if

            colat(k) = colat(k)*pi/180.d0
            elong(k) = elong(k)*pi/180.d0
            z(k) = h*1.d-3

        end do

        do k=1,10,1
            print *,colat(k),elong(k),z(k)
        end do
        
    end subroutine leer_archivo_dem

    subroutine procesar_datos_escribir_salida(pathout)

        implicit none

        character(len=100), intent(in) :: pathout
        
        open(unit=12,file=trim(pathout))
        
    end subroutine procesar_datos_escribir_salida
    
end module paleo_martesol_subs