program kinetics

! reaction A + B ---> C
!
! respective concentrations: a,b,c
! kinetic constant
! time step h
real*8 a,b,c,k

a=1d0
b=2d0
c=0d0
k=1d0
h=1d-2
t=0

open(unit=201,file='a.dat')
open(unit=202,file='b.dat')
open(unit=203,file='c.dat')


do i=1,100

    a=a-h*v(a,b)
    b=b-h*v(a,b)
    c=c+h*v(a,b)
    t=t+h

    write(201,*) t,a
    write(202,*) t,b
    write(203,*) t,c
    
enddo


contains

real function v(x,y)
real*8 x,y
v = k*x*y
end function

end program kinetics