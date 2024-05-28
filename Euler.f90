program Euler
Implicit none           ! liberdade para declarar todas as variáveis driblemando o ação automatica do programa de definir real ou integer dependendo da letra declarada 
real*8 :: a, b, c, d    ! abcd variáveis de concentração
real*8 :: dt, h, tr     ! dt variação no tempo, h passo, tr taxa de reação
integer :: n

a = 1.0
b = 1.0
h = 0.1
n = 10
g = 0.1

do i = 1, n
    c = c + tr * h ! preciso entender como incluir a participação de a e/ou b
    d = d + tr * h ! preciso entender como incluir a participação de a e/ou b
    a = a - tr * h ! preciso entender como a e b se relacionam
    b = b - tr * h ! preciso entender como a e b se relacionam
    write (*,*) a, b, c, d
End do
