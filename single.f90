program main
implicit none
real*8, dimension(5) :: v,w,u,k1,k2,k3,k4
real*8 tf,h,t,a
integer cont,cont2,i,j,k,cont3,crit,nrint,nrpoints,crit2,ni
integer conts,conts2
integer, parameter :: lk = selected_int_kind(18)
real*8 Cd,Ntot,Faraday,eappm,rm,cF,cFA,cWater,frt,ratek2
character(len=1) :: nfac
 Cd = 24d-6
ratek2=1.1d2
Ntot = 2.18d-9
Faraday = 9.648533d4
cF=3.55d-7

cFA = 1.99d-4
cWater = 55d-3
frt = Faraday/(8.314d0*298.15d0)

!!!rs=1d3
rm=5d2


h=1d-5
tf=1d2 !08h27

!nrint=1000000000_lk
nrint=int(t/h)
!write(6,*) 'nint',nrint
nrpoints=1000
crit2=2
crit=10000!int(nrint/nrpoints)
write(6,*) 'crit�rio',crit



do ni=1,1
write(6,*) 'j factor', ni

!!!kcpl = -9d2

!!!eapps0 = 9.1d-1
!eapps = eapps0
eappm = 8.9d-1

t=0d0
write(nfac,'(i1)') ni

!!!open(unit=201,file='es8'//nfac//'FB.dat')
!!!open(unit=202,file='es8'//nfac//'FL.dat')
!!!open(unit=203,file='es8'//nfac//'CO.dat') ! slave
!!!open(unit=204,file='es8'//nfac//'OH.dat')
!!!open(unit=205,file='es8'//nfac//'FI.dat')

open(unit=401,file='em9'//nfac//'FB.dat')
open(unit=402,file='em9'//nfac//'FL.dat')
open(unit=403,file='em9'//nfac//'CO.dat') ! master
open(unit=404,file='em9'//nfac//'OH.dat')
open(unit=405,file='em9'//nfac//'FI.dat')


cont=0
cont2=0

!v=(/5.84409d-2,1.45971d-4,5.79934d-1,3.28223d-4,5.80213d-1/)
v=(/0d0,0d0,0.58d0,0d0,0.05d0/)
!!!s=(/0d0,0d0,0.58d0,0d0,0.05d0/)

!v=(/2.1003413d-3,2.5981246d-4,6.55610424d-1,6.649863609d-6,2.64241242d-1/)

cont3=0
do k=1,4
cont3=cont3+1

do i=1,10000000!nrint
        call uf(v)                ! k1  --------------------------
        do j=1,size(u)
           k1(j)=u(j)             ! k1
        enddo
        
        do j=1,size(v)                  ! k2  ---------------------------
        w(j) = v(j) + h*k1(j)/2d0 
        enddo
        call uf(w)                ! k2
        do j=1,size(u)
           k2(j)=u(j)             ! k2
        enddo

        do j=1,size(v)                ! k3  ------------------------------
        w(j) = v(j) + h*k2(j)/2d0 
        enddo
        call uf(w)                ! k3
        do j=1,size(u)
           k3(j)=u(j)             ! k3
        enddo


        do j=1,size(v)                  ! k4   -----------------------------
        w(j) = v(j) + h*k3(j) 
        enddo
        call uf(w)                ! k4
        do j=1,size(u)
           k4(j)=u(j)             ! k4
        enddo

        do j=1,size(v)                ! new values
        v(j) = v(j) + h*(k1(j)+2d0*k2(j)+2d0*k3(j)+k4(j))/6d0 
        enddo 
       

        if (v(1).ne.v(1)) then
        write(11,*) cont3,i
        exit
        end if


        !write(6,*) cont!v(1),v(5)
        
        cont=cont+1
        t=t+h
        !if ((cont.eq.crit) .and. (cont3.eq.2)) then


        if (cont.eq.crit) then
                write(401,*) t,v(1)
                write(402,*) t,v(2)
                write(403,*) t,v(3)   !  MASTER
                write(404,*) t,v(4)
                write(405,*) t,v(5)
                cont2=cont2+1
                cont=0
                write(6,*) "n =",ni,"   saving  =",cont2

        end if


        ! master current calculation

        !!!eapps = eapps0 + kcpl * ( ((eappm-v(5))/rm) - ((eapps-s(5))/rm) )
        ! start slave integration


        !!!call ufs(s)                ! k1  --------------------------
        !!!do j=1,size(us)
        !!!     k1(j)=us(j)             ! k1
        !!!enddo

        !!!do j=1,size(s)                  ! k2  ---------------------------
        !!!w(j) = s(j) + h*k1(j)/2d0
        !!!enddo
        !!!call ufs(w)                ! k2
        !!!do j=1,size(us)
        !!!   k2(j)=us(j)             ! k2
        !!!enddo

        !!!do j=1,size(s)                ! k3  ------------------------------
        !!!w(j) = s(j) + h*k2(j)/2d0
        !!!enddo
        !!!call uf(w)                ! k3
        !!!do j=1,size(us)
        !!!   k3(j)=us(j)             ! k3
        !!!enddo


        !!!do j=1,size(s)                  ! k4   -----------------------------
        !!!w(j) = s(j) + h*k3(j)
        !!!enddo
        !!!call ufs(w)                ! k4
        !!!do j=1,size(us)
        !!!   k4(j)=us(j)             ! k4
        !!!enddo

        !!!do j=1,size(s)                ! new values
        !!!s(j) = s(j) + h*(k1(j)+2d0*k2(j)+2d0*k3(j)+k4(j))/6d0
        !!!enddo


        !!!if (s(1).ne.s(1)) then
        !!!write(12,*) cont3,i
        !!!exit
        !!!end if


        !write(6,*) cont!v(1),v(5)

        !!!conts=conts+1

        !if ((cont.eq.crit) .and. (cont3.eq.2)) then


        !!!if (conts.eq.crit) then
                !!!write(201,*) t,s(1)
                !!!write(202,*) t,s(2)
                !!!write(203,*) t,s(3)
                !!!write(204,*) t,s(4)
                !!!write(205,*) t,s(5)
                !!!conts2=conts2+1
                !!!conts=0
                !!!write(6,*) "n =",ni,"   saving  =",conts2

        !!!end if


        
        ! end slave integration
        


        
  enddo ! i ciclos de integração

enddo ! ni ciclo no valor da corrente

enddo

contains

       real function kox(beta,fi0,fi)
       real*8, intent(in):: beta,fi0,fi
       kox= dexp(beta*frt*(fi-fi0))
       end function

       real function kred(beta,fi0,fi)
       real*8, intent(in):: beta,fi0,fi
       kred= dexp(-(1d0-beta)*frt*(fi-fi0))
       end function

       real function v1(x)
               real*8, intent(in) :: x(:)
               v1 = cFA*vac(x)*vac(x)*kox(5d-1,-4d-2,x(5))
               end function


       real function v1r(x)
               real*8, intent(in) :: x(:)
               v1r = x(1)*kred(5d-1,2d-2,x(5))
               end function


       real function v2(x)
               real*8, intent(in) :: x(:)
               v2 = ratek2*x(1)*vac(x)*vac(x)
               end function

       real function v3(x)
               real*8, intent(in) :: x(:)
               v3 = x(2)*vac(x)*kox(5d-1,3.8d-1,x(5))
               end function

       real function v4(x)
               real*8, intent(in) :: x(:)
               v4 = x(2)*vac(x)*vac(x)*kred(5d-1,6d-1,x(5))
               end function


       real function v5(x)
               real*8, intent(in) :: x(:)
               v5 = cWater*vac(x)*kox(6.5d-1,4d-1,x(5))
               end function

       real function v5r(x)
               real*8, intent(in) :: x(:)
               v5r = x(4)*kred(6.5d-1,7.1d-1,x(5))
               end function


       real function v6(x)
               real*8, intent(in) :: x(:)
               v6 = x(3)*x(4)*kox(5d-1,7.8d-1,x(5))
               end function

       real function v7(x)
               real*8, intent(in) :: x(:)
               v7 = cF*vac(x)*kox(5d-1,-4d-1,x(5))
               end function


       real function f1(x)
               real*8, intent(in) :: x(:)
               f1 = v1(x) - v1r(x) - v2(x)
               end function

       real function f2(x)
               real*8, intent(in) :: x(:)
               f2 = v2(x) - v3(x) - v4(x)
               end function


       real function f3(x)
               real*8, intent(in) :: x(:)
               f3 = v4(x) - v6(x)
               end function

       real function f4(x)
               real*8, intent(in) :: x(:)
               f4 = v5(x) - v5r(x) - v6(x)
               end function

       real function f5(x)
               real*8, intent(in) :: x(:)
               f5 = ( (eappm-x(5))/rm - Faraday*Ntot*(v1(x)-v1r(x)+v3(x)-v4(x)+v5(x)-v5r(x)+v6(x)+2d0*v7(x)) )/Cd
               end function


       !!!real function f5s(x)
          !!!     real*8, intent(in) :: x(:)
             !!!  f5s = ( (eapps-x(5))/rs - Faraday*Ntot*(v1(x)-v1r(x)+v3(x)-v4(x)+v5(x)-v5r(x)+v6(x)+2d0*v7(x)) )/Cd
               !!!end function


       real function vac(x)
       real*8 x(:)
       vac=1d0-2d0*x(1)-x(2)-1.516d0*x(3)-x(4)
       end function
       
       !
       !  recFB x(1)
       !  recFL x(2)
       !  recCO x(3)
       !  recOH x(4)
       !
       
       
      subroutine uf(x)
      real*8 x(:)
      u=(/f1(x),f2(x),f3(x),f4(x),f5(x)/)
      end subroutine uf

      !!!subroutine ufs(x)
      !!!real*8 x(:)
      !!!us=(/f1(x),f2(x),f3(x),f4(x),f5s(x)/)
      !!!end subroutine ufs


end program
