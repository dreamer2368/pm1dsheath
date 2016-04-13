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

	subroutine twostream_initialize(this,Np,v0,vT,mode)		!generate initial distribution
		type(PM1D), intent(inout) :: this
		integer, intent(in) :: Np, mode
		real(mp), intent(in) :: v0, vT
		real(mp) :: xp0(Np), vp0(Np),rho_back
		real(mp) :: L,qs,ms
		integer :: i,j,N,pm(Np)
		L = this%L
		N = this%n

		qs = -this%wp*this%wp/(this%n*Np/L)
		ms = -qs
		rho_back = -qs*this%n*Np/L

		!velocity distribution initialize
		vp0 = vT*randn(Np)
		pm = (/ ( i, i=1,Np ) /)
		pm = 1 - 2*MOD(pm,2)
		vp0 = vp0 + pm*v0

		do i=1,this%n
			!spatial distribution initialize
			xp0 = (/ ( j*L/Np, j=0,Np-1 ) /) + 0.5*(i-1)*L/Np
			xp0 = xp0 + this%A0*L/Np*SIN( 2.0_mp*pi*xp0/L*mode )

			call buildSpecies(this%p(i),qs,ms)
			call setSpecies(this%p(i),Np,xp0,vp0)
		end do
		call setMesh(this%m,rho_back*(/ ( 1, i=1,this%m%ng) /))
	end subroutine

end module