program Euler
Implicit none          
real*8 :: A, B, t    
real*8 :: h, n, dt, A0, t0
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

Contains 

function 

End do
