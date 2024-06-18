module paleo_martesol_subs

    use,intrinsic :: iso_fortran_env, only : dp=>real64

    implicit none

    real(dp), parameter :: pi = 4.d0*datan(1.d0)
    real(dp), parameter :: radm = 3396.19d0
    real(dp), parameter :: uakm = 149597870.7d0
    real(dp), parameter :: epse = 23.43928d0*pi/180.d0

    real(dp) :: t0,dt,fdt
    
contains

    subroutine calcular_posicion_sol()

        implicit none

        character(len=100) :: filein,pathout

        real(dp), dimension(10000000) :: colat,elong,z

        ! Estrategia

        call leer_archivo_entrada(filein,pathout)

        call leer_archivo_dem(filein,colat,elong,z)

        call procesar_datos_escribir_salida()

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

    subroutine leer_archivo_dem(filein,colat,elong,z)
        
        implicit none

        integer :: k,feof,kmax

        real(dp) :: long,lat,h

        real(dp), dimension(:),intent(out) :: colat,elong,z

        character(len=100),intent(in) :: filein

        open(unit=11,file=trim(filein))

        k = 1

        
        feof = 0
        
        print *, 'procesando archivo "'//trim(filein)//'"'

        do while(feof.eq.0)

         read(11,*,iostat=feof) long,lat,h

         if (feof.gt.0) then
          print *,'revisar archivo de entrada'
          exit
         else if (feof.lt.0) then
          print *,  'listo!'
          exit
         else

          if (lat.ge.0.d0) then
           colat(k) = 90.d0 - lat
          else
           colat(k) = 90.d0 + dabs(lat)
          end if

          if (long.ge.0.d0) then
           elong(k) = long
          else
           elong(k) = 360.d0 - dabs(long)
          end if

          colat(k) = colat(k)*pi/180.d0
          elong(k) = elong(k)*pi/180.d0
              z(k) = h*1.d-3

         end if
         k = k + 1
        end do 
        kmax = k - 1
        
    end subroutine leer_archivo_dem

    subroutine procesar_datos_escribir_salida()

        implicit none

        
        
    end subroutine procesar_datos_escribir_salida
    
end module paleo_martesol_subs