module init

	use modPM1D
	use random
	implicit none

contains

	function randn(N) result(x)
		integer, intent(in) :: N
		real(mp) :: x(N)
		integer :: i

		do i = 1,N
			x(i) = random_normal()
		end do
	end function

	subroutine particle_initialize(this,v0,vT,mode,xp0,vp0,qs,ms,rho_back)		!generate initial distribution
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: mode
		real(mp), intent(in) :: v0, vT
		real(mp), intent(out) :: xp0(this%n), vp0(this%n), qs(this%n),ms(this%n),rho_back
		real(mp) :: L,qe,me
		integer :: i,N,pm(this%n)
		L = this%L
		N = this%n

		qe = -this%wp*this%wp/(this%n/L)
		qs = qe
		me = -qe
		ms = me
		rho_back = -qe*this%n/L

		!spatial distribution initialize
		xp0 = (/ ( i*L/N, i=0,N-1 ) /)
		xp0 = xp0 + this%A0*L/this%n*SIN( 2.0_mp*pi*xp0/L*mode )

		!velocity distribution initialize
		vp0 = vT*randn(N)
		pm = (/ ( i, i=1,N ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0 = vp0 + pm*v0
	end subroutine

end module