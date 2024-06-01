program kinetics

    ! reaction A + B ---> C
    !
    ! vector of concentrations: c[a,b,c]
    ! kinetic constant k
    ! time step h
    real*8 h,k
    real*8, dimension(3) :: c,r, re, ce ! adicionei re "r de Euler" e ce "C de euler"

    integer :: j
    real*8 :: t

    c=(/1d0,2d0,0d0/)
    k=1d0
    h=1d-2
    t=0d0

    open(unit=201,file='av.dat')
    open(unit=202,file='bv.dat')
    open(unit=203,file='cv.dat')

    ! Salva as concentraÃ§Ãµes iniciais
    write(201,*) t,c(1)
    write(202,*) t,c(2)
    write(203,*) t,c(3)

    do j=1,100

        call subf(c)
        re = r 
        ce = c + h * re  ! Estimativa inicial (M‚todo de Euler)
        
        call subf(ce)
        t = t + h
        do i=1,3
            c(i) = c(i) + h / 2 * (re(i) + r(i))
        end do

        ! Salva as concentraÃ§Ãµes atualizadas
        write(201,*) t,c(1)
        write(202,*) t,c(2)
        write(203,*) t,c(3)

    end do

    close(201)
    close(202)
    close(203)

contains

    real function vr(x)
        real*8 x(:)
        vr = -k*x(1)*x(2)
    end function


    real function vp(x)
        real*8 x(:)
        vp = k*x(1)*x(2)
    end function

    subroutine subf(x)
        real*8 x(:)
        r=(/vr(x),vr(x),vp(x)/)
    end subroutine subf

end program kinetics
