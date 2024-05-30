program Euler
    implicit none
    real(8) :: t, dt, A, B, C, k
    integer :: n, n_steps
    real(8), dimension(:), allocatable :: A_values, B_values, C_values
    character(len=100) :: output_filename
    integer :: iunit

    interface
        real(8) function taxa_reacao(conc_A, conc_B, k)
            real(8), intent(in) :: conc_A, conc_B, k
        end function taxa_reacao
    end interface

    ! Valores iniciais
    t = 0.0
    A = 1.0
    B = 2.0
    C = 0.0
    k = 0.1
    dt = 0.01
    n_steps = 100
    output_filename = "output_data.f90" ! Nome do arquivo de saÃ­da

    ! AlocaÃ§Ã£o de memÃ³ria para os arrays
    allocate(A_values(n_steps))
    allocate(B_values(n_steps))
    allocate(C_values(n_steps))

    ! AplicaÃ§Ã£o do mÃ©todo de Euler e armazenamento dos valores
    do n = 1, n_steps
        A_values(n) = A
        B_values(n) = B
        C_values(n) = C

        ! AtualizaÃ§Ã£o das concentraÃ§Ãµes de A e B
        A = A - dt * taxa_reacao(A, B, k)
        B = B - dt * taxa_reacao(A, B, k)
        C = C + dt * taxa_reacao(A, B, k)
        t = t + dt
    end do

    ! Salvar os dados em um arquivo
    open(unit=iunit, file=output_filename, status='replace')
    do n = 1, n_steps
        write(iunit, *) t, A_values(n), B_values(n), C_values(n)
    end do
    close(iunit)

    ! LiberaÃ§Ã£o da memÃ³ria alocada
    deallocate(A_values)
    deallocate(B_values)
    deallocate(C_values)

    stop
end program Euler

! FunÃ§Ã£o para a taxa de reaÃ§Ã£o
real(8) function taxa_reacao(conc_A, conc_B, k)
    real(8), intent(in) :: conc_A, conc_B, k
    taxa_reacao = k * conc_A * conc_B
end function taxa_reacao
