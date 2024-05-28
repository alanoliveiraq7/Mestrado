program Euler
Implicit none          
real*8 :: a, b, c, d   
real*8 :: dt, h, tr     
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
