program kinetics

! reaction A + B ---> C
!
! vector of concentrations: c[a,b,c]
! kinetic constant k
! time step h
real*8 h,k
real*8, dimension(3):: c,r

c=(/1d0,2d0,0d0/)
k=1d0
h=1d-2
t=0d0

open(unit=201,file='av.dat')
open(unit=202,file='bv.dat')
open(unit=203,file='cv.dat')


do i=1,3
    write(200+i,*) t,c(i)
enddo

do j=1,100

    call subf(c)
    t=t+h
    do i=1,3
        c(i)=c(i)+h*r(i)
        write(200+i,*) t,c(i)
    enddo

    
    
enddo


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