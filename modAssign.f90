module modAssign

	use modSpecies
	use modMesh

	implicit none

	type pmAssign
		integer :: np
		integer :: ng
		integer :: order

		integer, allocatable :: g(:,:)
		real(mp), allocatable :: frac(:,:)
		real(mp), allocatable :: h(:)
	end type

contains

	subroutine buildAssign(this,ng,order)
		type(pmAssign), intent(out) :: this
		integer, intent(in) :: ng, order

		this%ng = ng
		this%order = order
	end subroutine

	subroutine destroyAssign(this)
		type(pmAssign), intent(inout) :: this

		if( allocated(this%g) ) then
			deallocate(this%g)
		end if
		if( allocated(this%frac) ) then
			deallocate(this%frac)
		end if
		if( allocated(this%h) ) then
			deallocate(this%h)
		end if
	end subroutine

	subroutine assignMatrix(this,m,xp)
		type(pmAssign), intent(inout) :: this
		type(mesh), intent(inout) :: m
		real(mp), intent(inout) :: xp(:)

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
		real(mp), intent(inout) :: xp(:)
		integer :: i, np
		integer :: g1, gl, gr
		real(mp) :: fracl, fracr		!fraction for left grid point
		real(mp) :: h

		np = size(xp)
		this%np = np
		if( allocated(this%g) ) then
			deallocate(this%g)
		end if
		if( allocated(this%frac) ) then
			deallocate(this%frac)
		end if
		if( allocated(this%h) ) then
			deallocate(this%h)
		end if
		allocate(this%g(np,this%order+1))
		allocate(this%frac(np,this%order+1))
		allocate(this%h(np))

		!assignment matrix
		do i=1,this%np
			g1 = FLOOR(xp(i)/m%dx)+1
			gl = g1
			gr = gl+1

			h = xp(i)/m%dx - g1 + 1.0_mp
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
		real(mp), intent(inout) :: xp(:)
		integer :: i, np
		integer :: g		!nearest grid point
		integer :: gl		!left grid point
		integer :: gr		!right grid point
		real(mp) :: fracl		!fraction for left grid point
		real(mp) :: frac0		!fraction for nearest grid point
		real(mp) :: fracr		!fraction for right grid point
		real(mp) :: h			!distance from nearest grid point

		np = size(xp)
		this%np = np
		if( allocated(this%g) ) then
			deallocate(this%g)
		end if
		if( allocated(this%frac) ) then
			deallocate(this%frac)
		end if
		if( allocated(this%h) ) then
			deallocate(this%h)
		end if
		allocate(this%g(np,this%order+1))
		allocate(this%frac(np,this%order+1))
		allocate(this%h(np))
		!assignment matrix
		do i=1,this%np
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
		type(pmAssign), intent(inout) :: this(:)
		type(species), intent(inout) :: p(:)
		type(mesh), intent(inout) :: m
		integer :: i,ip

		m%rho = 0.0_mp
		do ip = 1, size(p)
			do i=1,p(ip)%np
				m%rho( this(ip)%g(i,:) ) = m%rho( this(ip)%g(i,:) ) + p(ip)%spwt*p(ip)%qs/m%dx*this(ip)%frac(i,:)
			end do
		end do
		m%rho = m%rho
	end subroutine

	subroutine forceAssign(this,p,m)
		type(pmAssign), intent(inout) :: this
		type(species), intent(inout) :: p
		type(mesh), intent(inout) :: m
		integer :: i

		p%Ep = 0.0_mp
		do i=1,this%np
			p%Ep(i) = sum( this%frac(i,:)*m%E(this%g(i,:)) )
		end do
	end subroutine

end module