module init

	use declaration
	use random
	implicit none

contains

	subroutine setup()
		!reference perturbation amplitude
		A = A0
		B = B0
		!time step parameter
		dt = 0.2_mp
		Nt = CEILING(T/dt)
		dt = T/Nt
		Ni = FLOOR(Ti/dt) + 1
		print *, 'Ni=',Ni,', Nt=',Nt,', dt=',dt

		!Grid setup
		dx = L/Ng
		xg = (/ ( i*L/Ng, i=0,Ng-1 ) /)
	end subroutine

	function randn(N) result(x)
		integer, intent(in) :: N
		real(mp) :: x(N)
		integer :: i

		do i = 1,N
			x(i) = SQRT(-2.0_mp*LOG(RAND()))*COS(2.0_mp*pi*RAND())
		end do
	end function

	subroutine particle_initialize()
		!spatial distribution initialize
		xp0 = (/ ( i*L/N, i=0,N-1 ) /)
!		xp0 = xp0 + A*L/N*SIN( 2.0_mp*pi*xp0/L*mode )
		xp0 = xp0 - A*L/N*0.3_mp*( SIN( 2.0_mp*pi*xp0/L*mode )	&
									+ SIN( 4.0_mp*pi*xp0/L*mode )	&
									+ SIN( 6.0_mp*pi*xp0/L*mode ) )

		!velocity distribution initialize
!		vp0 = vT*randn(N)
		vp0 = 0.0_mp
		pm = (/ ( i, i=1,N ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0 = vp0 + pm*v0
	end subroutine

end module