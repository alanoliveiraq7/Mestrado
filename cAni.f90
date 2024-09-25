program main
implicit none
real*8, dimension(6) :: v,w,u,k1,k2,k3,k4 ! Alterado para incluir mais uma variável (cAni)
real*8 tf,h,t,a,cAni
integer cont,cont2,i,j,k,cont3,crit,nrint,nrpoints,crit2,ni
integer conts,conts2
integer, parameter :: lk = selected_int_kind(18)
real*8 Cd,Ntot,Faraday,eappm,rm,cF,cFA,cWater,frt,ratek2,kcpl
character(len=1) :: nfac

Cd = 24d-6
ratek2 = 1.1d2
Ntot = 2.18d-9
Faraday = 9.648533d4
cF = 3.55d-7
cFA = 1.99d-4
cWater = 55d-3
frt = Faraday / (8.314d0 * 298.15d0)

rm = 2d3
cAni = 1d-5 ! Valor inicial para a nova espécie

h = 1d-5
tf = 1d2
nrint = int(t / h)
nrpoints = 1000
crit2 = 2
crit = 10000
write(6,*) 'critério', crit

do ni = 1, 1
    write(6,*) 'j factor', ni

    kcpl = -9d2
    eappm = 8.9d-1

    t = 0d0
    write(nfac,'(i1)') ni

    open(unit = 401, file = 'em9'//nfac//'FB.dat')
    open(unit = 402, file = 'em9'//nfac//'FL.dat')
    open(unit = 403, file = 'em9'//nfac//'CO.dat')
    open(unit = 404, file = 'em9'//nfac//'OH.dat')
    open(unit = 405, file = 'em9'//nfac//'FI.dat')
    open(unit = 406, file = 'em9'//nfac//'ANI.dat') ! Novo arquivo para cAni

    cont = 0
    cont2 = 0

    v = (/0d0, 0d0, 0.58d0, 0d0, 0.05d0, cAni/) ! Incluindo cAni no vetor de estados

    cont3 = 0
    do k = 1, 4
        cont3 = cont3 + 1

        do i = 1, 1000000
            call uf(v)
            do j = 1, size(u)
                k1(j) = u(j)
            enddo

            do j = 1, size(v)
                w(j) = v(j) + h * k1(j) / 2d0
            enddo
            call uf(w)
            do j = 1, size(u)
                k2(j) = u(j)
            enddo

            do j = 1, size(v)
                w(j) = v(j) + h * k2(j) / 2d0
            enddo
            call uf(w)
            do j = 1, size(u)
                k3(j) = u(j)
            enddo

            do j = 1, size(v)
                w(j) = v(j) + h * k3(j)
            enddo
            call uf(w)
            do j = 1, size(u)
                k4(j) = u(j)
            enddo

            do j = 1, size(v)
                v(j) = v(j) + h * (k1(j) + 2d0 * k2(j) + 2d0 * k3(j) + k4(j)) / 6d0
            enddo 

            if (v(1).ne.v(1)) then
                write(11,*) cont3,i
                exit
            end if

            cont = cont + 1
            t = t + h

            if (cont.eq.crit) then
                write(401,*) t,v(1)
                write(402,*) t,v(2)
                write(403,*) t,v(3)
                write(404,*) t,v(4)
                write(405,*) t,v(5)
                write(406,*) t,v(6) ! Salvando cAni no novo arquivo
                cont2 = cont2 + 1
                cont = 0
                write(6,*) "n =",ni,"   saving  =",cont2
            end if
        enddo
    enddo
enddo

contains

! Funções diretas e reversíveis para cAni
real function v8(x)
    real*8, intent(in) :: x(:)
    v8 = cAni * vac(x) * kox(4.5d-1, 3d-1, x(5)) ! Exemplo de comportamento direto
end function

real function v8r(x)
    real*8, intent(in) :: x(:)
    v8r = x(6) * kred(4.5d-1, 6d-1, x(5)) ! Exemplo de comportamento reversível
end function

! Atualização de uf para incluir a nova função
subroutine uf(x)
    real*8 x(:)
    u = (/f1(x), f2(x), f3(x), f4(x), f5(x), v8(x) - v8r(x)/) ! Atualizando as equações para incluir cAni
end subroutine uf

! Outras funções continuam sem alterações

end program