module modAssign

	use modPlasma
	use modMesh

	implicit none

	type pmAssign
		integer :: n
		integer :: ng
		integer :: order

		integer, allocatable :: g(:,:)
		real(mp), allocatable :: frac(:,:)
		real(mp), allocatable :: h(:)
	end type

contains

	subroutine buildAssign(this,n,ng,order)
		type(pmAssign), intent(out) :: this
		integer, intent(in) :: n, ng, order

		this%n = n
		this%ng = ng
		this%order = order

		allocate(this%g(n,order+1))
		allocate(this%frac(n,order+1))
		allocate(this%h(n))
	end subroutine

	subroutine destroyAssign(this)
		type(pmAssign), intent(inout) :: this

		deallocate(this%g)
		deallocate(this%frac)
		deallocate(this%h)
	end subroutine

	subroutine assignMatrix(this,m,xp)
		type(pmAssign), intent(inout) :: this
		type(mesh), intent(inout) :: m
		real(mp), intent(inout) :: xp(this%n)

		SELECT CASE(this%order)
			CASE(1)
				call assign_CIC(this,m,xp)
			CASE(2)
				call assign_TSC(this,m,xp)
		END SELECT
	end subroutine

	subroutine assign_CIC(this,m,xp)	!apply BC and create assignment matrix
		type(pmAssign), intent(inout) :: this
		type(mesh), intent(inout) :: m
		real(mp), intent(inout) :: xp(this%n)
		integer :: i
		integer :: g1, gl, gr
		real(mp) :: fracl, fracr		!fraction for left grid point
		real(mp) :: h

		!assignment matrix
		do i=1,this%n
			g1 = FLOOR(xp(i)/m%dx - 0.5_mp)+1
			gl = g1
			gr = gl+1

			h = xp(i)/m%dx - g1 + 0.5_mp
			fracl = 1.0_mp - ABS(h)
			fracr = 1.0_mp - fracl

			this%h(i) = h
			this%g(i,:) = (/ gl, gr /)
			this%frac(i,:) = (/ fracl, fracr /)
		end do
	end subroutine

	subroutine assign_TSC(this,m,xp)	!apply BC and create assignment matrix
		type(pmAssign), intent(inout) :: this
		type(mesh), intent(inout) :: m
		real(mp), intent(inout) :: xp(this%n)
		integer :: i
		integer :: g		!nearest grid point
		integer :: gl		!left grid point
		integer :: gr		!right grid point
		real(mp) :: fracl		!fraction for left grid point
		real(mp) :: frac0		!fraction for nearest grid point
		real(mp) :: fracr		!fraction for right grid point
		real(mp) :: h			!distance from nearest grid point

		!assignment matrix
		do i=1,this%n
			g = FLOOR(xp(i)/m%dx + 0.5_mp)+1
			gl = g-1
			gr = g+1
			h = xp(i)/m%dx - g + 1.0_mp

			frac0 = 0.75_mp - h*h
			fracl = 0.5_mp*(0.5_mp-h)*(0.5_mp-h)
			fracr = 0.5_mp*(0.5_mp+h)*(0.5_mp+h)

			this%g(i,:) = (/gl,g,gr/)
			this%h(i) = h
			this%frac(i,:) = (/ fracl, frac0, fracr /)
		end do
	end subroutine

	subroutine chargeAssign(this,p,m)
		type(pmAssign), intent(inout) :: this
		type(plasma), intent(inout) :: p
		type(mesh), intent(inout) :: m
		integer :: i

		m%rho = 0.0_mp
		do i=1,this%n
			m%rho( this%g(i,:) ) = m%rho( this%g(i,:) ) + p%qs(i)/m%dx*this%frac(i,:)
		end do
		m%rho = m%rho + m%rho_back
	end subroutine

	subroutine forceAssign(this,p,m)
		type(pmAssign), intent(inout) :: this
		type(plasma), intent(inout) :: p
		type(mesh), intent(inout) :: m
		integer :: i

		p%Ep = 0.0_mp
		do i=1,this%n
			p%Ep(i) = sum( this%frac(i,:)*m%E(this%g(i,:)) )
		end do
	end subroutine

end module